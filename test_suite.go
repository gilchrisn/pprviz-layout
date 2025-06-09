// test_suite.go - Comprehensive test suite for GraphProcessor

package main

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"reflect"
	"sort"
	"strings"
	"testing"
	"time"
	"github.com/gilchrisn/pprviz-layout/graphprocessor"
)

// TestResult tracks test outcomes
type TestResult struct {
	Name     string
	Passed   bool
	Error    error
	Duration time.Duration
	Details  string
}

// TestSuite manages all tests
type TestSuite struct {
	Results    []TestResult
	TempDir    string
	Processor  *graphprocessor.GraphProcessor
	PassCount  int
	FailCount  int
}

// NewTestSuite creates a new test suite
func NewTestSuite() *TestSuite {
	tempDir := filepath.Join(os.TempDir(), fmt.Sprintf("graphprocessor_test_%d", time.Now().Unix()))
	os.MkdirAll(tempDir, 0755)
	
	config := graphprocessor.DefaultConfig()
	config.Verbose = false // Quiet during tests
	
	return &TestSuite{
		TempDir:   tempDir,
		Processor: graphprocessor.NewGraphProcessor(tempDir, config),
		Results:   make([]TestResult, 0),
	}
}

// RunTest executes a single test and records results
func (ts *TestSuite) RunTest(name string, testFunc func() error) {
	start := time.Now()
	err := testFunc()
	duration := time.Since(start)
	
	result := TestResult{
		Name:     name,
		Passed:   err == nil,
		Error:    err,
		Duration: duration,
	}
	
	if err == nil {
		ts.PassCount++
		fmt.Printf("✓ %s (%.2fms)\n", name, float64(duration.Nanoseconds())/1e6)
	} else {
		ts.FailCount++
		fmt.Printf("✗ %s: %v (%.2fms)\n", name, err, float64(duration.Nanoseconds())/1e6)
	}
	
	ts.Results = append(ts.Results, result)
}

// Cleanup removes temporary files
func (ts *TestSuite) Cleanup() {
	os.RemoveAll(ts.TempDir)
}

// PrintSummary prints test results summary
func (ts *TestSuite) PrintSummary() {
	fmt.Printf("\n=== Test Summary ===\n")
	fmt.Printf("Total: %d, Passed: %d, Failed: %d\n", 
		len(ts.Results), ts.PassCount, ts.FailCount)
	
	if ts.FailCount > 0 {
		fmt.Printf("\nFailed Tests:\n")
		for _, result := range ts.Results {
			if !result.Passed {
				fmt.Printf("  - %s: %v\n", result.Name, result.Error)
			}
		}
	}
}

// Test data generators
func (ts *TestSuite) createTestGraph(name string, edges [][]int, n int) error {
	datasetDir := filepath.Join(ts.TempDir, "dataset")
	if err := os.MkdirAll(datasetDir, 0755); err != nil {
		return err
	}

	// Create graph file
	graphFile := filepath.Join(datasetDir, name+".txt")
	content := ""
	for _, edge := range edges {
		if len(edge) == 2 {
			content += fmt.Sprintf("%d %d\n", edge[0], edge[1])
		}
	}
	
	if err := os.WriteFile(graphFile, []byte(content), 0644); err != nil {
		return err
	}

	// Create attribute file
	attrFile := filepath.Join(datasetDir, name+"_attribute.txt")
	attrContent := fmt.Sprintf("n=%d\nm=%d\n", n, len(edges))
	return os.WriteFile(attrFile, []byte(attrContent), 0644)
}

func (ts *TestSuite) createTestHierarchy(name string, k int, mapping map[string][]int, hierarchy map[string][]int, roots []string) error {
	mappingDir := filepath.Join(ts.TempDir, "louvain", "mapping-output")
	hierarchyDir := filepath.Join(ts.TempDir, "louvain", "hierachy-output")
	
	if err := os.MkdirAll(mappingDir, 0755); err != nil {
		return err
	}
	if err := os.MkdirAll(hierarchyDir, 0755); err != nil {
		return err
	}

	// Create mapping file
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

	// Create hierarchy file
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
	rootContent := strings.Join(roots, "\n") + "\n"
	return os.WriteFile(rootFile, []byte(rootContent), 0644)
}

// ==== GRAPH LOADING TESTS ====

func (ts *TestSuite) testEmptyGraph() error {
	if err := ts.createTestGraph("empty", [][]int{}, 0); err != nil {
		return err
	}
	
	err := ts.Processor.LoadGraph("empty")
	if err == nil {
		return fmt.Errorf("expected error for empty graph")
	}
	return nil // Error expected
}

func (ts *TestSuite) testSingleNodeGraph() error {
	if err := ts.createTestGraph("single", [][]int{}, 1); err != nil {
		return err
	}
	
	return ts.Processor.LoadGraph("single")
}

func (ts *TestSuite) testSelfLoopGraph() error {
	edges := [][]int{{0, 0}, {0, 1}, {1, 2}}
	if err := ts.createTestGraph("selfloop", edges, 3); err != nil {
		return err
	}
	
	if err := ts.Processor.LoadGraph("selfloop"); err != nil {
		return err
	}
	
	// Verify self-loops are ignored
	stats := ts.Processor.GetGraphStats()
	if edges, ok := stats["edges"].(int64); !ok || edges != 2 {
		return fmt.Errorf("expected 2 edges (self-loop ignored), got %v", edges)
	}
	
	return nil
}

func (ts *TestSuite) testDisconnectedGraph() error {
	edges := [][]int{{0, 1}, {2, 3}} // Two disconnected components
	if err := ts.createTestGraph("disconnected", edges, 4); err != nil {
		return err
	}
	
	return ts.Processor.LoadGraph("disconnected")
}

func (ts *TestSuite) testMalformedGraphFile() error {
	datasetDir := filepath.Join(ts.TempDir, "dataset")
	os.MkdirAll(datasetDir, 0755)
	
	graphFile := filepath.Join(datasetDir, "malformed.txt")
	malformedContent := "0 1\ninvalid line\n2\n3 4 5\n"
	os.WriteFile(graphFile, []byte(malformedContent), 0644)
	
	attrFile := filepath.Join(datasetDir, "malformed_attribute.txt")
	os.WriteFile(attrFile, []byte("n=5\nm=2\n"), 0644)
	
	return ts.Processor.LoadGraph("malformed") // Should handle gracefully
}

func (ts *TestSuite) testTriangleGraph() error {
	edges := [][]int{{0, 1}, {1, 2}, {2, 0}}
	if err := ts.createTestGraph("triangle", edges, 3); err != nil {
		return err
	}
	
	if err := ts.Processor.LoadGraph("triangle"); err != nil {
		return err
	}
	
	// Verify graph properties
	stats := ts.Processor.GetGraphStats()
	if nodes, ok := stats["nodes"].(int); !ok || nodes != 3 {
		return fmt.Errorf("expected 3 nodes, got %v", nodes)
	}
	if edges, ok := stats["edges"].(int64); !ok || edges != 3 {
		return fmt.Errorf("expected 3 edges, got %v", edges)
	}
	
	return nil
}

// ==== HIERARCHY LOADING TESTS ====

func (ts *TestSuite) testSimpleHierarchy() error {
	mapping := map[string][]int{
		"c0_l1_0": {0, 1},
		"c0_l1_1": {2},
	}
	hierarchy := map[string][]int{
		"c0_l2_0": {0, 1},
	}
	roots := []string{"c0_l2_0"}
	
	if err := ts.createTestHierarchy("simple", 25, mapping, hierarchy, roots); err != nil {
		return err
	}
	
	return ts.Processor.LoadHierarchy("simple", 25)
}

func (ts *TestSuite) testEmptyHierarchy() error {
	mapping := map[string][]int{}
	hierarchy := map[string][]int{}
	roots := []string{}
	
	if err := ts.createTestHierarchy("emptyhier", 25, mapping, hierarchy, roots); err != nil {
		return err
	}
	
	return ts.Processor.LoadHierarchy("emptyhier", 25)
}

func (ts *TestSuite) testInconsistentHierarchy() error {
	mapping := map[string][]int{
		"c0_l1_0": {0, 1},
		"c0_l1_1": {2, 999}, // Node 999 doesn't exist
	}
	hierarchy := map[string][]int{
		"c0_l2_0": {0, 1},
	}
	roots := []string{"c0_l2_0"}
	
	if err := ts.createTestHierarchy("inconsistent", 25, mapping, hierarchy, roots); err != nil {
		return err
	}
	
	return ts.Processor.LoadHierarchy("inconsistent", 25) // Should handle gracefully
}

// ==== DNPR BUILDING TESTS ====

func (ts *TestSuite) testDNPRTriangle() error {
	edges := [][]int{{0, 1}, {1, 2}, {2, 0}}
	if err := ts.createTestGraph("dnpr_triangle", edges, 3); err != nil {
		return err
	}
	
	if err := ts.Processor.LoadGraph("dnpr_triangle"); err != nil {
		return err
	}
	
	if err := ts.Processor.BuildDNPRIndices("dnpr_triangle"); err != nil {
		return err
	}
	
	// Verify DNPR file was created
	dnprPath := filepath.Join(ts.TempDir, "pr_idx", "dnpr_triangle.dnpr")
	if _, err := os.Stat(dnprPath); os.IsNotExist(err) {
		return fmt.Errorf("DNPR file was not created")
	}
	
	return nil
}

func (ts *TestSuite) testDNPRStar() error {
	// Star graph: center node connected to all others
	edges := [][]int{{0, 1}, {0, 2}, {0, 3}, {0, 4}}
	if err := ts.createTestGraph("dnpr_star", edges, 5); err != nil {
		return err
	}
	
	if err := ts.Processor.LoadGraph("dnpr_star"); err != nil {
		return err
	}
	
	return ts.Processor.BuildDNPRIndices("dnpr_star")
}

func (ts *TestSuite) testDNPRPath() error {
	// Path graph: 0-1-2-3-4
	edges := [][]int{{0, 1}, {1, 2}, {2, 3}, {3, 4}}
	if err := ts.createTestGraph("dnpr_path", edges, 5); err != nil {
		return err
	}
	
	if err := ts.Processor.LoadGraph("dnpr_path"); err != nil {
		return err
	}
	
	return ts.Processor.BuildDNPRIndices("dnpr_path")
}

func (ts *TestSuite) testDNPRSingleNode() error {
	if err := ts.createTestGraph("dnpr_single", [][]int{}, 1); err != nil {
		return err
	}
	
	if err := ts.Processor.LoadGraph("dnpr_single"); err != nil {
		return err
	}
	
	return ts.Processor.BuildDNPRIndices("dnpr_single")
}

// ==== PDIST BUILDING TESTS ====

func (ts *TestSuite) testPDistWithoutDNPR() error {
	edges := [][]int{{0, 1}, {1, 2}, {2, 0}}
	if err := ts.createTestGraph("pdist_no_dnpr", edges, 3); err != nil {
		return err
	}
	
	mapping := map[string][]int{"c0_l1_0": {0, 1, 2}}
	hierarchy := map[string][]int{"c0_l2_0": {0}}
	roots := []string{"c0_l2_0"}
	
	if err := ts.createTestHierarchy("pdist_no_dnpr", 25, mapping, hierarchy, roots); err != nil {
		return err
	}
	
	if err := ts.Processor.LoadGraph("pdist_no_dnpr"); err != nil {
		return err
	}
	
	if err := ts.Processor.LoadHierarchy("pdist_no_dnpr", 25); err != nil {
		return err
	}
	
	// Should fail without DNPR
	err := ts.Processor.BuildPDistIndices("pdist_no_dnpr", 25)
	if err == nil {
		return fmt.Errorf("expected error when building PDist without DNPR")
	}
	
	return nil // Error expected
}

func (ts *TestSuite) testPDistComplete() error {
	edges := [][]int{{0, 1}, {1, 2}, {2, 0}}
	if err := ts.createTestGraph("pdist_complete", edges, 3); err != nil {
		return err
	}
	
	mapping := map[string][]int{"c0_l1_0": {0, 1, 2}}
	hierarchy := map[string][]int{"c0_l2_0": {0}}
	roots := []string{"c0_l2_0"}
	
	if err := ts.createTestHierarchy("pdist_complete", 25, mapping, hierarchy, roots); err != nil {
		return err
	}
	
	if err := ts.Processor.LoadGraph("pdist_complete"); err != nil {
		return err
	}
	
	if err := ts.Processor.LoadHierarchy("pdist_complete", 25); err != nil {
		return err
	}
	
	if err := ts.Processor.BuildDNPRIndices("pdist_complete"); err != nil {
		return err
	}
	
	return ts.Processor.BuildPDistIndices("pdist_complete", 25)
}

// ==== COORDINATE GENERATION TESTS ====

func (ts *TestSuite) testCoordinatesLevel1() error {
	edges := [][]int{{0, 1}, {1, 2}, {2, 3}}
	if err := ts.createTestGraph("coord_l1", edges, 4); err != nil {
		return err
	}
	
	mapping := map[string][]int{
		"c0_l1_0": {0, 1, 2, 3},
	}
	hierarchy := map[string][]int{
		"c0_l2_0": {0},
	}
	roots := []string{"c0_l2_0"}
	
	if err := ts.createTestHierarchy("coord_l1", 25, mapping, hierarchy, roots); err != nil {
		return err
	}
	
	if err := ts.Processor.LoadGraph("coord_l1"); err != nil {
		return err
	}
	
	if err := ts.Processor.LoadHierarchy("coord_l1", 25); err != nil {
		return err
	}
	
	if err := ts.Processor.BuildDNPRIndices("coord_l1"); err != nil {
		return err
	}
	
	if err := ts.Processor.BuildPDistIndices("coord_l1", 25); err != nil {
		return err
	}
	
	coords, err := ts.Processor.GetSupernodeCoordinates("coord_l1", "c0_l1_0", 25)
	if err != nil {
		return err
	}
	
	// Validate coordinate structure
	if len(coords.X) != 4 || len(coords.Y) != 4 || len(coords.Radii) != 4 {
		return fmt.Errorf("expected 4 coordinates, got X:%d Y:%d R:%d", 
			len(coords.X), len(coords.Y), len(coords.Radii))
	}
	
	if coords.Metadata.Level != 1 {
		return fmt.Errorf("expected level 1, got %d", coords.Metadata.Level)
	}
	
	return nil
}

func (ts *TestSuite) testCoordinatesInvalidSupernode() error {
	edges := [][]int{{0, 1}}
	if err := ts.createTestGraph("coord_invalid", edges, 2); err != nil {
		return err
	}
	
	mapping := map[string][]int{"c0_l1_0": {0, 1}}
	hierarchy := map[string][]int{"c0_l2_0": {0}}
	roots := []string{"c0_l2_0"}
	
	if err := ts.createTestHierarchy("coord_invalid", 25, mapping, hierarchy, roots); err != nil {
		return err
	}
	
	if err := ts.Processor.LoadGraph("coord_invalid"); err != nil {
		return err
	}
	
	if err := ts.Processor.LoadHierarchy("coord_invalid", 25); err != nil {
		return err
	}
	
	if err := ts.Processor.BuildDNPRIndices("coord_invalid"); err != nil {
		return err
	}
	
	if err := ts.Processor.BuildPDistIndices("coord_invalid", 25); err != nil {
		return err
	}
	
	_, err := ts.Processor.GetSupernodeCoordinates("coord_invalid", "nonexistent", 25)
	if err == nil {
		return fmt.Errorf("expected error for invalid supernode")
	}
	
	return nil // Error expected
}

func (ts *TestSuite) testCoordinatesSinglePoint() error {
	if err := ts.createTestGraph("coord_single", [][]int{}, 1); err != nil {
		return err
	}
	
	mapping := map[string][]int{"c0_l1_0": {0}}
	hierarchy := map[string][]int{"c0_l2_0": {0}}
	roots := []string{"c0_l2_0"}
	
	if err := ts.createTestHierarchy("coord_single", 25, mapping, hierarchy, roots); err != nil {
		return err
	}
	
	if err := ts.Processor.LoadGraph("coord_single"); err != nil {
		return err
	}
	
	if err := ts.Processor.LoadHierarchy("coord_single", 25); err != nil {
		return err
	}
	
	if err := ts.Processor.BuildDNPRIndices("coord_single"); err != nil {
		return err
	}
	
	if err := ts.Processor.BuildPDistIndices("coord_single", 25); err != nil {
		return err
	}
	
	coords, err := ts.Processor.GetSupernodeCoordinates("coord_single", "c0_l1_0", 25)
	if err != nil {
		return err
	}
	
	if len(coords.X) != 1 {
		return fmt.Errorf("expected 1 coordinate for single node, got %d", len(coords.X))
	}
	
	return nil
}

// ==== MDS ALGORITHM TESTS ====

func (ts *TestSuite) testMDSIdenticalPoints() error {
	// Test MDS with identical distance matrix (should handle gracefully)
	processor := graphprocessor.NewGraphProcessor(ts.TempDir, graphprocessor.DefaultConfig())
	
	// Create a simple 2x2 distance matrix with identical points
	distMatrix := [][]float64{
		{0, 0},
		{0, 0},
	}
	
	// This should not crash
	coords, err := processor.PerformMDS(distMatrix)
	if err != nil {
		return err
	}
	
	if len(coords) != 2 || len(coords[0]) != 2 || len(coords[1]) != 2 {
		return fmt.Errorf("unexpected coordinate dimensions")
	}
	
	return nil
}

func (ts *TestSuite) testMDSValidTriangle() error {
	// Test MDS with a valid triangle distance matrix
	processor := graphprocessor.NewGraphProcessor(ts.TempDir, graphprocessor.DefaultConfig())
	
	// Create a 3x3 distance matrix for an equilateral triangle
	distMatrix := [][]float64{
		{0, 1, 1},
		{1, 0, 1},
		{1, 1, 0},
	}
	
	coords, err := processor.performMDS(distMatrix)
	if err != nil {
		return err
	}
	
	if len(coords) != 2 || len(coords[0]) != 3 || len(coords[1]) != 3 {
		return fmt.Errorf("unexpected coordinate dimensions")
	}
	
	// Verify distances are approximately preserved
	for i := 0; i < 3; i++ {
		for j := i + 1; j < 3; j++ {
			dx := coords[0][i] - coords[0][j]
			dy := coords[1][i] - coords[1][j]
			actualDist := math.Sqrt(dx*dx + dy*dy)
			expectedDist := distMatrix[i][j]
			
			if math.Abs(actualDist-expectedDist) > 0.1 {
				return fmt.Errorf("distance not preserved: expected %.2f, got %.2f", 
					expectedDist, actualDist)
			}
		}
	}
	
	return nil
}

// ==== FILE I/O TESTS ====

func (ts *TestSuite) testFileIOCorruption() error {
	edges := [][]int{{0, 1}}
	if err := ts.createTestGraph("fileio", edges, 2); err != nil {
		return err
	}
	
	if err := ts.Processor.LoadGraph("fileio"); err != nil {
		return err
	}
	
	if err := ts.Processor.BuildDNPRIndices("fileio"); err != nil {
		return err
	}
	
	// Corrupt the DNPR file
	dnprPath := filepath.Join(ts.TempDir, "pr_idx", "fileio.dnpr")
	if err := os.WriteFile(dnprPath, []byte("corrupted"), 0644); err != nil {
		return err
	}
	
	// Try to load corrupted file
	_, err := ts.Processor.LoadDNPR(dnprPath)
	if err == nil {
		return fmt.Errorf("expected error when loading corrupted DNPR file")
	}
	
	return nil // Error expected
}

func (ts *TestSuite) testSaveLoadCoordinates() error {
	coords := &graphprocessor.CoordinateResult{
		X:     []float64{1.0, 2.0, 3.0},
		Y:     []float64{4.0, 5.0, 6.0},
		Radii: []float64{0.1, 0.2, 0.3},
		Metadata: &graphprocessor.SupernodeMetadata{
			SupernodeName: "test",
			Level:         1,
			Children:      []string{"a", "b", "c"},
		},
	}
	
	filename := filepath.Join(ts.TempDir, "test_coords.json")
	
	if err := ts.Processor.SaveCoordinates(coords, filename); err != nil {
		return err
	}
	
	loaded, err := ts.Processor.LoadCoordinates(filename)
	if err != nil {
		return err
	}
	
	if !reflect.DeepEqual(coords.X, loaded.X) ||
		!reflect.DeepEqual(coords.Y, loaded.Y) ||
		!reflect.DeepEqual(coords.Radii, loaded.Radii) {
		return fmt.Errorf("coordinates not preserved after save/load")
	}
	
	return nil
}

// ==== STRESS TESTS ====

func (ts *TestSuite) testLargeGraph() error {
	// Create a larger graph for stress testing
	n := 100
	edges := make([][]int, 0)
	
	// Create a ring graph
	for i := 0; i < n; i++ {
		edges = append(edges, []int{i, (i + 1) % n})
	}
	
	if err := ts.createTestGraph("large", edges, n); err != nil {
		return err
	}
	
	start := time.Now()
	if err := ts.Processor.LoadGraph("large"); err != nil {
		return err
	}
	loadTime := time.Since(start)
	
	if loadTime > 5*time.Second {
		return fmt.Errorf("graph loading took too long: %v", loadTime)
	}
	
	start = time.Now()
	if err := ts.Processor.BuildDNPRIndices("large"); err != nil {
		return err
	}
	buildTime := time.Since(start)
	
	if buildTime > 10*time.Second {
		return fmt.Errorf("DNPR building took too long: %v", buildTime)
	}
	
	return nil
}

func (ts *TestSuite) testConcurrentAccess() error {
	edges := [][]int{{0, 1}, {1, 2}, {2, 0}}
	if err := ts.createTestGraph("concurrent", edges, 3); err != nil {
		return err
	}
	
	if err := ts.Processor.LoadGraph("concurrent"); err != nil {
		return err
	}
	
	// Test concurrent stats access
	done := make(chan bool, 10)
	for i := 0; i < 10; i++ {
		go func() {
			stats := ts.Processor.GetGraphStats()
			if stats["nodes"] != 3 {
				panic("concurrent access failed")
			}
			done <- true
		}()
	}
	
	for i := 0; i < 10; i++ {
		<-done
	}
	
	return nil
}

// Main test runner
func main() {
	ts := NewTestSuite()
	defer ts.Cleanup()
	
	fmt.Println("Starting Graph Processor Test Suite...")
	fmt.Println("=====================================")
	
	// Graph Loading Tests
	fmt.Println("\n--- Graph Loading Tests ---")
	ts.RunTest("Empty Graph", ts.testEmptyGraph)
	ts.RunTest("Single Node Graph", ts.testSingleNodeGraph)
	ts.RunTest("Self Loop Graph", ts.testSelfLoopGraph)
	ts.RunTest("Disconnected Graph", ts.testDisconnectedGraph)
	ts.RunTest("Malformed Graph File", ts.testMalformedGraphFile)
	ts.RunTest("Triangle Graph", ts.testTriangleGraph)
	
	// Hierarchy Loading Tests
	fmt.Println("\n--- Hierarchy Loading Tests ---")
	ts.RunTest("Simple Hierarchy", ts.testSimpleHierarchy)
	ts.RunTest("Empty Hierarchy", ts.testEmptyHierarchy)
	ts.RunTest("Inconsistent Hierarchy", ts.testInconsistentHierarchy)
	
	// DNPR Building Tests
	fmt.Println("\n--- DNPR Building Tests ---")
	ts.RunTest("DNPR Triangle", ts.testDNPRTriangle)
	ts.RunTest("DNPR Star", ts.testDNPRStar)
	ts.RunTest("DNPR Path", ts.testDNPRPath)
	ts.RunTest("DNPR Single Node", ts.testDNPRSingleNode)
	
	// PDist Building Tests
	fmt.Println("\n--- PDist Building Tests ---")
	ts.RunTest("PDist Without DNPR", ts.testPDistWithoutDNPR)
	ts.RunTest("PDist Complete", ts.testPDistComplete)
	
	// Coordinate Generation Tests
	fmt.Println("\n--- Coordinate Generation Tests ---")
	ts.RunTest("Coordinates Level 1", ts.testCoordinatesLevel1)
	ts.RunTest("Coordinates Invalid Supernode", ts.testCoordinatesInvalidSupernode)
	ts.RunTest("Coordinates Single Point", ts.testCoordinatesSinglePoint)
	
	// MDS Algorithm Tests
	fmt.Println("\n--- MDS Algorithm Tests ---")
	ts.RunTest("MDS Identical Points", ts.testMDSIdenticalPoints)
	ts.RunTest("MDS Valid Triangle", ts.testMDSValidTriangle)
	
	// File I/O Tests
	fmt.Println("\n--- File I/O Tests ---")
	ts.RunTest("File I/O Corruption", ts.testFileIOCorruption)
	ts.RunTest("Save/Load Coordinates", ts.testSaveLoadCoordinates)
	
	// Stress Tests
	fmt.Println("\n--- Stress Tests ---")
	ts.RunTest("Large Graph", ts.testLargeGraph)
	ts.RunTest("Concurrent Access", ts.testConcurrentAccess)
	
	ts.PrintSummary()
	
	if ts.FailCount > 0 {
		os.Exit(1)
	}
}