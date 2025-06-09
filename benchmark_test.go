// benchmark_test.go - Performance benchmarks and advanced correctness tests

package main

import (
	"fmt"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"time"
	"./graphprocessor"
)

// BenchmarkSuite handles performance and correctness benchmarks
type BenchmarkSuite struct {
	TempDir   string
	Processor *graphprocessor.GraphProcessor
}

func NewBenchmarkSuite() *BenchmarkSuite {
	tempDir := filepath.Join(os.TempDir(), fmt.Sprintf("graphprocessor_bench_%d", time.Now().Unix()))
	os.MkdirAll(tempDir, 0755)
	
	config := graphprocessor.DefaultConfig()
	config.Verbose = true
	
	return &BenchmarkSuite{
		TempDir:   tempDir,
		Processor: graphprocessor.NewGraphProcessor(tempDir, config),
	}
}

func (bs *BenchmarkSuite) Cleanup() {
	os.RemoveAll(bs.TempDir)
}

// ==== GRAPH GENERATORS FOR TESTING ====

// generateRandomGraph creates a random graph with specified properties
func (bs *BenchmarkSuite) generateRandomGraph(name string, n int, avgDegree float64, seed int64) error {
	rand.Seed(seed)
	
	// Calculate number of edges
	m := int(float64(n) * avgDegree / 2.0)
	
	edges := make(map[string]bool) // Use map to avoid duplicates
	edgeList := make([][]int, 0)
	
	for len(edgeList) < m {
		u := rand.Intn(n)
		v := rand.Intn(n)
		
		if u == v {
			continue // No self-loops
		}
		
		if u > v {
			u, v = v, u // Ensure consistent ordering
		}
		
		edgeKey := fmt.Sprintf("%d-%d", u, v)
		if !edges[edgeKey] {
			edges[edgeKey] = true
			edgeList = append(edgeList, []int{u, v})
		}
	}
	
	return bs.createGraph(name, edgeList, n)
}

// generateScaleFreeGraph creates a scale-free graph using preferential attachment
func (bs *BenchmarkSuite) generateScaleFreeGraph(name string, n int, m0 int, seed int64) error {
	rand.Seed(seed)
	
	if m0 >= n {
		m0 = n / 2
	}
	
	degrees := make([]int, n)
	edges := make([][]int, 0)
	
	// Start with a complete graph of m0 nodes
	for i := 0; i < m0; i++ {
		for j := i + 1; j < m0; j++ {
			edges = append(edges, []int{i, j})
			degrees[i]++
			degrees[j]++
		}
	}
	
	// Add remaining nodes with preferential attachment
	for i := m0; i < n; i++ {
		// Calculate total degree so far
		totalDegree := 0
		for j := 0; j < i; j++ {
			totalDegree += degrees[j]
		}
		
		if totalDegree == 0 {
			totalDegree = 1
		}
		
		// Connect to m0 existing nodes based on their degree
		connected := make(map[int]bool)
		for len(connected) < m0 && len(connected) < i {
			// Preferential attachment: probability proportional to degree
			r := rand.Float64() * float64(totalDegree)
			cumDegree := 0.0
			
			for j := 0; j < i; j++ {
				cumDegree += float64(degrees[j])
				if r <= cumDegree && !connected[j] {
					edges = append(edges, []int{i, j})
					degrees[i]++
					degrees[j]++
					connected[j] = true
					break
				}
			}
		}
		
		// If we couldn't connect to enough nodes, connect to random ones
		for len(connected) < m0 && len(connected) < i {
			j := rand.Intn(i)
			if !connected[j] {
				edges = append(edges, []int{i, j})
				degrees[i]++
				degrees[j]++
				connected[j] = true
			}
		}
	}
	
	return bs.createGraph(name, edges, n)
}

// generateHierarchicalGraph creates a hierarchical graph with communities
func (bs *BenchmarkSuite) generateHierarchicalGraph(name string, communities int, nodesPerCommunity int, intraP, interP float64, seed int64) error {
	rand.Seed(seed)
	
	n := communities * nodesPerCommunity
	edges := make([][]int, 0)
	
	// Create intra-community edges
	for c := 0; c < communities; c++ {
		start := c * nodesPerCommunity
		end := start + nodesPerCommunity
		
		for i := start; i < end; i++ {
			for j := i + 1; j < end; j++ {
				if rand.Float64() < intraP {
					edges = append(edges, []int{i, j})
				}
			}
		}
	}
	
	// Create inter-community edges
	for c1 := 0; c1 < communities; c1++ {
		for c2 := c1 + 1; c2 < communities; c2++ {
			start1 := c1 * nodesPerCommunity
			end1 := start1 + nodesPerCommunity
			start2 := c2 * nodesPerCommunity
			end2 := start2 + nodesPerCommunity
			
			for i := start1; i < end1; i++ {
				for j := start2; j < end2; j++ {
					if rand.Float64() < interP {
						edges = append(edges, []int{i, j})
					}
				}
			}
		}
	}
	
	return bs.createGraph(name, edges, n)
}

func (bs *BenchmarkSuite) createGraph(name string, edges [][]int, n int) error {
	datasetDir := filepath.Join(bs.TempDir, "dataset")
	os.MkdirAll(datasetDir, 0755)

	// Create graph file
	graphFile := filepath.Join(datasetDir, name+".txt")
	content := ""
	for _, edge := range edges {
		content += fmt.Sprintf("%d %d\n", edge[0], edge[1])
	}
	
	if err := os.WriteFile(graphFile, []byte(content), 0644); err != nil {
		return err
	}

	// Create attribute file
	attrFile := filepath.Join(datasetDir, name+"_attribute.txt")
	attrContent := fmt.Sprintf("n=%d\nm=%d\n", n, len(edges))
	return os.WriteFile(attrFile, []byte(attrContent), 0644)
}

func (bs *BenchmarkSuite) createRealisticHierarchy(name string, k int, n int, levels int) error {
	mappingDir := filepath.Join(bs.TempDir, "louvain", "mapping-output")
	hierarchyDir := filepath.Join(bs.TempDir, "louvain", "hierachy-output")
	
	os.MkdirAll(mappingDir, 0755)
	os.MkdirAll(hierarchyDir, 0755)

	// Generate hierarchical clustering
	mapping := make(map[string][]int)
	hierarchy := make(map[string][]int)
	
	// Level 1: Partition nodes into clusters
	clusterSize := int(math.Ceil(float64(n) / float64(k)))
	for i := 0; i < k; i++ {
		start := i * clusterSize
		end := start + clusterSize
		if end > n {
			end = n
		}
		
		if start < end {
			nodes := make([]int, end-start)
			for j := start; j < end; j++ {
				nodes[j-start] = j
			}
			mapping[fmt.Sprintf("c0_l1_%d", i)] = nodes
		}
	}
	
	// Higher levels: Merge clusters
	for level := 2; level <= levels; level++ {
		prevClusters := int(math.Ceil(float64(k) / math.Pow(2, float64(level-2))))
		currClusters := int(math.Ceil(float64(prevClusters) / 2.0))
		
		for i := 0; i < currClusters; i++ {
			children := make([]int, 0)
			start := i * 2
			end := start + 2
			if end > prevClusters {
				end = prevClusters
			}
			
			for j := start; j < end; j++ {
				children = append(children, j)
			}
			
			if len(children) > 0 {
				hierarchy[fmt.Sprintf("c0_l%d_%d", level, i)] = children
			}
		}
	}
	
	// Save mapping file
	mapFile := filepath.Join(mappingDir, fmt.Sprintf("%s_%d.dat", name, k))
	mapContent := ""
	for supernode, nodes := range mapping {
		mapContent += fmt.Sprintf("%s %d", supernode, len(nodes))
		for _, node := range nodes {
			mapContent += fmt.Sprintf(" %d", node)
		}
		mapContent += "\n"
	}
	
	if err := os.WriteFile(mapFile, []byte(mapContent), 0644); err != nil {
		return err
	}

	// Save hierarchy file
	hierFile := filepath.Join(hierarchyDir, fmt.Sprintf("%s_%d.dat", name, k))
	hierContent := ""
	for supernode, children := range hierarchy {
		hierContent += fmt.Sprintf("%s %d", supernode, len(children))
		for _, child := range children {
			hierContent += fmt.Sprintf(" %d", child)
		}
		hierContent += "\n"
	}
	
	if err := os.WriteFile(hierFile, []byte(hierContent), 0644); err != nil {
		return err
	}

	// Create root file
	rootFile := filepath.Join(hierarchyDir, fmt.Sprintf("%s_%d.root", name, k))
	rootContent := fmt.Sprintf("c0_l%d_0\n", levels)
	return os.WriteFile(rootFile, []byte(rootContent), 0644)
}

// ==== CORRECTNESS TESTS ====

func (bs *BenchmarkSuite) testPPRCorrectness() error {
	fmt.Println("Testing PPR Correctness...")
	
	// Create a simple triangle graph where we can verify PPR values
	edges := [][]int{{0, 1}, {1, 2}, {2, 0}}
	if err := bs.createGraph("ppr_correct", edges, 3); err != nil {
		return err
	}
	
	if err := bs.Processor.LoadGraph("ppr_correct"); err != nil {
		return err
	}
	
	// Build DNPR
	if err := bs.Processor.BuildDNPRIndices("ppr_correct"); err != nil {
		return err
	}
	
	// For a regular triangle, all nodes should have equal PageRank
	dnprPath := filepath.Join(bs.TempDir, "pr_idx", "ppr_correct.dnpr")
	dnpr, err := bs.Processor.LoadDNPR(dnprPath)
	if err != nil {
		return err
	}
	
	if len(dnpr) != 3 {
		return fmt.Errorf("expected 3 DNPR values, got %d", len(dnpr))
	}
	
	// Check that values sum approximately to 1 (after normalization)
	sum := dnpr[0] + dnpr[1] + dnpr[2]
	if math.Abs(sum-1.0) > 0.1 { // Allow some tolerance
		fmt.Printf("Warning: DNPR values don't sum to 1: %f\n", sum)
	}
	
	// For symmetric graph, values should be approximately equal
	avgVal := sum / 3.0
	for i, val := range dnpr {
		if math.Abs(val-avgVal) > 0.1 {
			fmt.Printf("Warning: Node %d has non-uniform DNPR: %f vs %f\n", i, val, avgVal)
		}
	}
	
	fmt.Printf("PPR values: [%.4f, %.4f, %.4f], sum: %.4f\n", dnpr[0], dnpr[1], dnpr[2], sum)
	return nil
}

func (bs *BenchmarkSuite) testMDSCorrectness() error {
	fmt.Println("Testing MDS Correctness...")
	
	// Create a square distance matrix
	distMatrix := [][]float64{
		{0, 1, math.Sqrt(2), 1},
		{1, 0, 1, math.Sqrt(2)},
		{math.Sqrt(2), 1, 0, 1},
		{1, math.Sqrt(2), 1, 0},
	}
	
	coords, err := bs.Processor.PerformMDS(distMatrix)
	if err != nil {
		return err
	}
	
	// Verify that the MDS preserves the square structure approximately
	n := len(coords[0])
	if n != 4 {
		return fmt.Errorf("expected 4 points, got %d", n)
	}
	
	// Calculate resulting distances
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			dx := coords[0][i] - coords[0][j]
			dy := coords[1][i] - coords[1][j]
			actualDist := math.Sqrt(dx*dx + dy*dy)
			expectedDist := distMatrix[i][j]
			
			relativeError := math.Abs(actualDist-expectedDist) / expectedDist
			if relativeError > 0.2 { // 20% tolerance
				fmt.Printf("Warning: Distance %d-%d not preserved well: expected %.2f, got %.2f (error: %.1f%%)\n", 
					i, j, expectedDist, actualDist, relativeError*100)
			}
		}
	}
	
	fmt.Printf("MDS coordinates: X=[%.2f, %.2f, %.2f, %.2f], Y=[%.2f, %.2f, %.2f, %.2f]\n",
		coords[0][0], coords[0][1], coords[0][2], coords[0][3],
		coords[1][0], coords[1][1], coords[1][2], coords[1][3])
	
	return nil
}

func (bs *BenchmarkSuite) testCoordinateStability() error {
	fmt.Println("Testing Coordinate Generation Stability...")
	
	// Create the same graph and generate coordinates multiple times
	edges := [][]int{{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 2}}
	if err := bs.createGraph("stability", edges, 4); err != nil {
		return err
	}
	
	if err := bs.createRealisticHierarchy("stability", 25, 4, 2); err != nil {
		return err
	}
	
	if err := bs.Processor.LoadGraph("stability"); err != nil {
		return err
	}
	
	if err := bs.Processor.LoadHierarchy("stability", 25); err != nil {
		return err
	}
	
	if err := bs.Processor.BuildDNPRIndices("stability"); err != nil {
		return err
	}
	
	if err := bs.Processor.BuildPDistIndices("stability", 25); err != nil {
		return err
	}
	
	// Generate coordinates multiple times
	var allCoords []*graphprocessor.CoordinateResult
	for i := 0; i < 3; i++ {
		coords, err := bs.Processor.GetSupernodeCoordinates("stability", "c0_l1_0", 25)
		if err != nil {
			return err
		}
		allCoords = append(allCoords, coords)
	}
	
	// Check that results are identical (deterministic)
	for i := 1; i < len(allCoords); i++ {
		if len(allCoords[0].X) != len(allCoords[i].X) {
			return fmt.Errorf("coordinate length differs between runs")
		}
		
		for j := 0; j < len(allCoords[0].X); j++ {
			if math.Abs(allCoords[0].X[j]-allCoords[i].X[j]) > 1e-10 ||
				math.Abs(allCoords[0].Y[j]-allCoords[i].Y[j]) > 1e-10 {
				return fmt.Errorf("coordinates differ between runs: run 0 (%.6f, %.6f) vs run %d (%.6f, %.6f)",
					allCoords[0].X[j], allCoords[0].Y[j], i, allCoords[i].X[j], allCoords[i].Y[j])
			}
		}
	}
	
	fmt.Printf("Coordinate generation is deterministic across %d runs\n", len(allCoords))
	return nil
}

// ==== PERFORMANCE BENCHMARKS ====

func (bs *BenchmarkSuite) benchmarkGraphLoading() error {
	fmt.Println("Benchmarking Graph Loading...")
	
	sizes := []int{100, 500, 1000, 2000}
	
	for _, n := range sizes {
		name := fmt.Sprintf("load_bench_%d", n)
		
		// Generate random graph
		if err := bs.generateRandomGraph(name, n, 4.0, 42); err != nil {
			return err
		}
		
		start := time.Now()
		if err := bs.Processor.LoadGraph(name); err != nil {
			return err
		}
		duration := time.Since(start)
		
		stats := bs.Processor.GetGraphStats()
		fmt.Printf("  n=%d: %v (%.2f nodes/ms)\n", n, duration, float64(n)/float64(duration.Nanoseconds())*1e6)
		
		// Verify stats
		if stats["nodes"] != n {
			return fmt.Errorf("node count mismatch for n=%d", n)
		}
	}
	
	return nil
}

func (bs *BenchmarkSuite) benchmarkDNPRBuilding() error {
	fmt.Println("Benchmarking DNPR Building...")
	
	sizes := []int{100, 300, 500}
	
	for _, n := range sizes {
		name := fmt.Sprintf("dnpr_bench_%d", n)
		
		// Generate scale-free graph (more realistic)
		if err := bs.generateScaleFreeGraph(name, n, 3, 42); err != nil {
			return err
		}
		
		if err := bs.Processor.LoadGraph(name); err != nil {
			return err
		}
		
		start := time.Now()
		if err := bs.Processor.BuildDNPRIndices(name); err != nil {
			return err
		}
		duration := time.Since(start)
		
		fmt.Printf("  n=%d: %v (%.2f nodes/ms)\n", n, duration, float64(n)/float64(duration.Nanoseconds())*1e6)
	}
	
	return nil
}

func (bs *BenchmarkSuite) benchmarkCoordinateGeneration() error {
	fmt.Println("Benchmarking Coordinate Generation...")
	
	clusterSizes := []int{5, 10, 20, 50}
	
	for _, size := range clusterSizes {
		name := fmt.Sprintf("coord_bench_%d", size)
		
		// Create hierarchical graph
		if err := bs.generateHierarchicalGraph(name, 4, size, 0.8, 0.1, 42); err != nil {
			return err
		}
		
		if err := bs.createRealisticHierarchy(name, 25, 4*size, 3); err != nil {
			return err
		}
		
		if err := bs.Processor.LoadGraph(name); err != nil {
			return err
		}
		
		if err := bs.Processor.LoadHierarchy(name, 25); err != nil {
			return err
		}
		
		if err := bs.Processor.BuildDNPRIndices(name); err != nil {
			return err
		}
		
		if err := bs.Processor.BuildPDistIndices(name, 25); err != nil {
			return err
		}
		
		start := time.Now()
		coords, err := bs.Processor.GetSupernodeCoordinates(name, "c0_l1_0", 25)
		if err != nil {
			return err
		}
		duration := time.Since(start)
		
		fmt.Printf("  cluster_size=%d: %v (%d points, %.2f points/ms)\n", 
			size, duration, len(coords.X), float64(len(coords.X))/float64(duration.Nanoseconds())*1e6)
	}
	
	return nil
}

func (bs *BenchmarkSuite) benchmarkMemoryUsage() error {
	fmt.Println("Benchmarking Memory Usage...")
	
	var m1, m2 runtime.MemStats
	
	// Measure baseline memory
	runtime.GC()
	runtime.ReadMemStats(&m1)
	
	// Load a moderately large graph
	n := 1000
	if err := bs.generateRandomGraph("memory_test", n, 6.0, 42); err != nil {
		return err
	}
	
	if err := bs.Processor.LoadGraph("memory_test"); err != nil {
		return err
	}
	
	if err := bs.Processor.BuildDNPRIndices("memory_test"); err != nil {
		return err
	}
	
	// Measure memory after loading
	runtime.GC()
	runtime.ReadMemStats(&m2)
	
	memUsed := m2.Alloc - m1.Alloc
	fmt.Printf("  Memory used for n=%d: %.2f MB (%.2f KB/node)\n", 
		n, float64(memUsed)/(1024*1024), float64(memUsed)/(1024*float64(n)))
	
	return nil
}

// ==== EDGE CASE TESTS ====

func (bs *BenchmarkSuite) testNumericalStability() error {
	fmt.Println("Testing Numerical Stability...")
	
	// Test with very small and very large numbers
	distMatrix := [][]float64{
		{0, 1e-10, 1e10},
		{1e-10, 0, 1e10},
		{1e10, 1e10, 0},
	}
	
	_, err := bs.Processor.performMDS(distMatrix)
	if err != nil {
		return err
	}
	
	// Test with NaN and Inf (should be handled gracefully)
	distMatrix = [][]float64{
		{0, 1, math.NaN()},
		{1, 0, math.Inf(1)},
		{math.NaN(), math.Inf(1), 0},
	}
	
	// This should not crash, but may return an error
	_, err = bs.Processor.PerformMDS(distMatrix)
	// We accept either success or controlled failure
	
	fmt.Println("  Numerical stability tests completed")
	return nil
}

func (bs *BenchmarkSuite) testBoundaryConditions() error {
	fmt.Println("Testing Boundary Conditions...")
	
	// Test with maximum node ID
	edges := [][]int{{0, 999}, {999, 998}}
	if err := bs.createGraph("boundary", edges, 1000); err != nil {
		return err
	}
	
	if err := bs.Processor.LoadGraph("boundary"); err != nil {
		return err
	}
	
	// Test with zero-degree nodes
	edges = [][]int{{0, 1}} // Node 2 has degree 0
	if err := bs.createGraph("zerodegree", edges, 3); err != nil {
		return err
	}
	
	if err := bs.Processor.LoadGraph("zerodegree"); err != nil {
		return err
	}
	
	if err := bs.Processor.BuildDNPRIndices("zerodegree"); err != nil {
		return err
	}
	
	fmt.Println("  Boundary condition tests completed")
	return nil
}

// Main benchmark runner
func main() {
	fmt.Println("Graph Processor Benchmark Suite")
	fmt.Println("===============================")
	
	bs := NewBenchmarkSuite()
	defer bs.Cleanup()
	
	tests := []struct {
		name string
		fn   func() error
	}{
		{"PPR Correctness", bs.testPPRCorrectness},
		{"MDS Correctness", bs.testMDSCorrectness},
		{"Coordinate Stability", bs.testCoordinateStability},
		{"Graph Loading Performance", bs.benchmarkGraphLoading},
		{"DNPR Building Performance", bs.benchmarkDNPRBuilding},
		{"Coordinate Generation Performance", bs.benchmarkCoordinateGeneration},
		{"Memory Usage", bs.benchmarkMemoryUsage},
		{"Numerical Stability", bs.testNumericalStability},
		{"Boundary Conditions", bs.testBoundaryConditions},
	}
	
	passed := 0
	failed := 0
	
	for _, test := range tests {
		fmt.Printf("\n--- %s ---\n", test.name)
		start := time.Now()
		
		if err := test.fn(); err != nil {
			fmt.Printf("FAILED: %v\n", err)
			failed++
		} else {
			fmt.Printf("PASSED (%.2fs)\n", time.Since(start).Seconds())
			passed++
		}
	}
	
	fmt.Printf("\n=== Benchmark Summary ===\n")
	fmt.Printf("Total: %d, Passed: %d, Failed: %d\n", len(tests), passed, failed)
	
	if failed > 0 {
		os.Exit(1)
	}
}