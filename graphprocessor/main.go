package graphprocessor

import (
	"bufio"
	// "encoding/gob"
	// "encoding/json"
	// "errors"
	"fmt"
	// "io"
	"math"
	// "math/rand"
	"os"
	"path/filepath"
	// "sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

// Config holds configuration parameters
type Config struct {
	Alpha     float64 // damping factor
	Delta     float64 // precision parameter
	PFail     float64 // failure probability
	EpR       float64 // epsilon parameter
	Tau       float64 // threshold parameter
	K         int     // cluster size
	ThreadNum int     // number of threads
	Verbose   bool    // verbose output
}

// DefaultConfig returns default configuration
func DefaultConfig() *Config {
	return &Config{
		Alpha:     0.2,
		Delta:     1.0 / (250.0 * math.Log(1000)), // will be updated based on graph size
		PFail:     0.001,                          // will be updated based on graph size
		EpR:       0.5,
		Tau:       1.0 / math.Sqrt(25*1000), // will be updated based on k and n
		K:         25,
		ThreadNum: 1,
		Verbose:   false,
	}
}

// Graph represents the graph structure
type Graph struct {
	Nodes   []int     // node list
	Edges   [][]int   // adjacency list
	Degrees []int     // node degrees
	N       int       // number of nodes
	M       int64     // number of edges
	DBar    float64   // average degree
	MaxLevel int      // maximum hierarchy level
	config  *Config
}

// SupernodeHierarchy holds the hierarchical clustering information
type SupernodeHierarchy struct {
	Super2Leaf    map[string][]int    // supernode to leaf nodes mapping
	Super2Super   map[string][]int    // supernode to child supernodes mapping
	Level1Cluster [][]int             // base level clusters
	Root          []string            // root supernodes
	HubCluster    []string            // hub clusters
}

// PPRIndex stores precomputed PPR indices
type PPRIndex struct {
	RandomWalkIdx       []int                    // random walk destinations
	RWIdxInfoOffset     []uint64                 // offsets for RW index
	RWIdxInfoSize       []uint64                 // sizes for RW index
	BackwardIdxTarget   []int                    // backward index targets
	BackwardIdxInCluster []float64              // backward index values
	BackwardIdxInfo     map[int][2]int          // target -> (offset, size)
	DegreePageRank      []float64               // degree-based PageRank
}

// ForwardIndex represents forward search results
type ForwardIndex struct {
	Reserve  map[int]float64 // reserve values
	Residual map[int]float64 // residual values
}

// CoordinateResult holds the coordinate generation results
type CoordinateResult struct {
	X        []float64             `json:"x"`
	Y        []float64             `json:"y"`
	Radii    []float64             `json:"radii"`
	Metadata *SupernodeMetadata    `json:"metadata"`
}

// SupernodeMetadata contains metadata about the supernode
type SupernodeMetadata struct {
	SupernodeName string     `json:"supernode_name"`
	Level         int        `json:"level"`
	Children      []string   `json:"children"`
	NodeWeights   []int      `json:"node_weights"`
	Degrees       []float64  `json:"degrees"`
	DPRValues     []float64  `json:"dpr_values"`
	LeafNodes     []int      `json:"leaf_nodes,omitempty"`
	PPRMatrix     [][]float64 `json:"ppr_matrix"`
	PDistMatrix   [][]float64 `json:"pdist_matrix"`
}

// GraphProcessor is the main processor interface
type GraphProcessor struct {
	config     *Config
	basePath   string
	graph      *Graph
	hierarchy  *SupernodeHierarchy
	pprIndex   *PPRIndex
	mu         sync.RWMutex
}

// NewGraphProcessor creates a new graph processor
func NewGraphProcessor(basePath string, config *Config) *GraphProcessor {
	if config == nil {
		config = DefaultConfig()
	}
	return &GraphProcessor{
		config:   config,
		basePath: basePath,
	}
}

// LoadGraph loads a graph from the dataset
func (gp *GraphProcessor) LoadGraph(datasetID string) error {
	gp.mu.Lock()
	defer gp.mu.Unlock()

	// Load graph structure
	graphPath := filepath.Join(gp.basePath, "dataset", datasetID+".txt")
	attrPath := filepath.Join(gp.basePath, "dataset", datasetID+"_attribute.txt")

	// Read attributes first
	n, m, err := gp.readAttributes(attrPath)
	if err != nil {
		return fmt.Errorf("failed to read attributes: %v", err)
	}

	// Initialize graph
	graph := &Graph{
		N:       n,
		M:       m,
		Nodes:   make([]int, n),
		Edges:   make([][]int, n),
		Degrees: make([]int, n),
		config:  gp.config,
	}

	// Initialize nodes
	for i := 0; i < n; i++ {
		graph.Nodes[i] = i
	}

	// Read edges
	if err := gp.readEdges(graphPath, graph); err != nil {
		return fmt.Errorf("failed to read edges: %v", err)
	}

	// Calculate average degree
	graph.DBar = float64(2*graph.M) / float64(graph.N)

	// Update config parameters based on graph size
	gp.config.Delta = 1.0 / (250.0 * math.Log(float64(graph.N)))
	gp.config.PFail = 1.0 / float64(graph.N)
	gp.config.Tau = 1.0 / math.Sqrt(float64(gp.config.K)*float64(graph.N))

	gp.graph = graph

	if gp.config.Verbose {
		fmt.Printf("Loaded graph: %d nodes, %d edges, avg degree: %.2f\n", 
			graph.N, graph.M, graph.DBar)
	}

	return nil
}

// LoadHierarchy loads the hierarchical clustering information
func (gp *GraphProcessor) LoadHierarchy(datasetID string, k int) error {
	gp.mu.Lock()
	defer gp.mu.Unlock()

	hierarchy := &SupernodeHierarchy{
		Super2Leaf:  make(map[string][]int),
		Super2Super: make(map[string][]int),
		Root:        make([]string, 0),
	}

	// Load mapping file (super2leaf)
	mapPath := filepath.Join(gp.basePath, "louvain", "mapping-output", 
		fmt.Sprintf("%s_%d.dat", datasetID, k))
	if err := gp.loadMapping(mapPath, hierarchy); err != nil {
		return fmt.Errorf("failed to load mapping: %v", err)
	}

	// Load hierarchy file (super2super)
	hierPath := filepath.Join(gp.basePath, "louvain", "hierachy-output", 
		fmt.Sprintf("%s_%d.dat", datasetID, k))
	maxLevel, err := gp.loadHierarchyFile(hierPath, hierarchy)
	if err != nil {
		return fmt.Errorf("failed to load hierarchy: %v", err)
	}

	// Load root file
	rootPath := filepath.Join(gp.basePath, "louvain", "hierachy-output", 
		fmt.Sprintf("%s_%d.root", datasetID, k))
	if err := gp.loadRoot(rootPath, hierarchy); err != nil {
		return fmt.Errorf("failed to load root: %v", err)
	}

	gp.hierarchy = hierarchy
	if gp.graph != nil {
		gp.graph.MaxLevel = maxLevel
	}

	if gp.config.Verbose {
		fmt.Printf("Loaded hierarchy: %d supernodes, %d level1 clusters, max level: %d\n", 
			len(hierarchy.Super2Leaf), len(hierarchy.Level1Cluster), maxLevel)
	}

	return nil
}

// readAttributes reads the graph attributes file
func (gp *GraphProcessor) readAttributes(filename string) (int, int64, error) {
	file, err := os.Open(filename)
	if err != nil {
		return 0, 0, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	var n int
	var m int64

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if strings.Contains(line, "=") {
			parts := strings.Split(line, "=")
			if len(parts) == 2 {
				key := strings.TrimSpace(parts[0])
				value := strings.TrimSpace(parts[1])
				
				if key == "n" {
					n, err = strconv.Atoi(value)
					if err != nil {
						return 0, 0, err
					}
				} else if key == "m" {
					m, err = strconv.ParseInt(value, 10, 64)
					if err != nil {
						return 0, 0, err
					}
				}
			}
		}
	}

	return n, m, scanner.Err()
}

// readEdges reads the graph edges file
func (gp *GraphProcessor) readEdges(filename string, graph *Graph) error {
	file, err := os.Open(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		parts := strings.Fields(line)
		if len(parts) != 2 {
			continue
		}

		u, err1 := strconv.Atoi(parts[0])
		v, err2 := strconv.Atoi(parts[1])

		if err1 != nil || err2 != nil {
			continue
		}

		if u >= graph.N || v >= graph.N || u == v {
			continue
		}

		// Add edges (undirected)
		graph.Edges[u] = append(graph.Edges[u], v)
		graph.Edges[v] = append(graph.Edges[v], u)
	}

	// Update degrees
	for i := 0; i < graph.N; i++ {
		graph.Degrees[i] = len(graph.Edges[i])
	}

	return scanner.Err()
}

// loadMapping loads the supernode to leaf mapping
func (gp *GraphProcessor) loadMapping(filename string, hierarchy *SupernodeHierarchy) error {
	file, err := os.Open(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		parts := strings.Fields(line)
		if len(parts) < 2 {
			continue
		}

		name := parts[0]
		size, err := strconv.Atoi(parts[1])
		if err != nil {
			continue
		}

		if len(parts) < 2+size {
			continue
		}

		nodes := make([]int, size)
		for i := 0; i < size; i++ {
			nodes[i], err = strconv.Atoi(parts[2+i])
			if err != nil {
				return err
			}
		}

		hierarchy.Super2Leaf[name] = nodes

		// Extract level information
		_, level := gp.extractSupernodeInfo(name)
		if level == 1 {
			hierarchy.Level1Cluster = append(hierarchy.Level1Cluster, nodes)
		}
	}

	return scanner.Err()
}

// loadHierarchyFile loads the supernode hierarchy file
func (gp *GraphProcessor) loadHierarchyFile(filename string, hierarchy *SupernodeHierarchy) (int, error) {
	file, err := os.Open(filename)
	if err != nil {
		return 0, err
	}
	defer file.Close()

	maxLevel := 0
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		parts := strings.Fields(line)
		if len(parts) < 2 {
			continue
		}

		name := parts[0]
		size, err := strconv.Atoi(parts[1])
		if err != nil {
			continue
		}

		if len(parts) < 2+size {
			continue
		}

		children := make([]int, size)
		for i := 0; i < size; i++ {
			children[i], err = strconv.Atoi(parts[2+i])
			if err != nil {
				return 0, err
			}
		}

		hierarchy.Super2Super[name] = children

		// Track maximum level
		_, level := gp.extractSupernodeInfo(name)
		if level > maxLevel {
			maxLevel = level
		}
	}

	return maxLevel, scanner.Err()
}

// loadRoot loads the root supernodes
func (gp *GraphProcessor) loadRoot(filename string, hierarchy *SupernodeHierarchy) error {
	file, err := os.Open(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line != "" {
			hierarchy.Root = append(hierarchy.Root, line)
		}
	}

	return scanner.Err()
}

// extractSupernodeInfo extracts component and level from supernode name
// Format: c{component}_l{level}_{id}
func (gp *GraphProcessor) extractSupernodeInfo(supernode string) (int, int) {
	// Parse supernode name like "c0_l2_5"
	parts := strings.Split(supernode, "_")
	if len(parts) != 3 {
		return 0, 0
	}

	component := 0
	level := 0

	if strings.HasPrefix(parts[0], "c") {
		if c, err := strconv.Atoi(parts[0][1:]); err == nil {
			component = c
		}
	}

	if strings.HasPrefix(parts[1], "l") {
		if l, err := strconv.Atoi(parts[1][1:]); err == nil {
			level = l
		}
	}

	return component, level
}

// BuildDNPRIndices builds the degree-normalized PageRank indices
func (gp *GraphProcessor) BuildDNPRIndices(datasetID string) error {
	if gp.graph == nil {
		if err := gp.LoadGraph(datasetID); err != nil {
			return err
		}
	}

	if gp.config.Verbose {
		fmt.Println("Building DNPR indices...")
	}

	start := time.Now()

	// Build degree-normalized PageRank
	dnpr, err := gp.buildDNPR(20) // 20 iterations
	if err != nil {
		return fmt.Errorf("failed to build DNPR: %v", err)
	}

	// Save DNPR indices
	dnprPath := filepath.Join(gp.basePath, "pr_idx", datasetID+".dnpr")
	if err := gp.saveDNPR(dnprPath, dnpr); err != nil {
		return fmt.Errorf("failed to save DNPR: %v", err)
	}

	if gp.config.Verbose {
		fmt.Printf("DNPR indices built in %v\n", time.Since(start))
	}

	return nil
}

// BuildPDistIndices builds the personalized distance indices  
func (gp *GraphProcessor) BuildPDistIndices(datasetID string, k int) error {
	if gp.graph == nil {
		if err := gp.LoadGraph(datasetID); err != nil {
			return err
		}
	}

	if gp.hierarchy == nil {
		if err := gp.LoadHierarchy(datasetID, k); err != nil {
			return err
		}
	}

	if gp.config.Verbose {
		fmt.Println("Building PDist indices...")
	}

	start := time.Now()

	// Load DNPR if not already loaded
	if gp.pprIndex == nil || len(gp.pprIndex.DegreePageRank) == 0 {
		dnprPath := filepath.Join(gp.basePath, "pr_idx", datasetID+".dnpr")
		dnpr, err := gp.loadDNPR(dnprPath)
		if err != nil {
			return fmt.Errorf("failed to load DNPR: %v", err)
		}
		if gp.pprIndex == nil {
			gp.pprIndex = &PPRIndex{}
		}
		gp.pprIndex.DegreePageRank = dnpr
	}

	// Build backward push indices
	if err := gp.buildBackwardPushIndex(); err != nil {
		return fmt.Errorf("failed to build backward push index: %v", err)
	}

	// Save backward indices
	bwdPath := filepath.Join(gp.basePath, "bwd_idx", datasetID)
	if err := gp.saveBackwardIndex(bwdPath); err != nil {
		return fmt.Errorf("failed to save backward index: %v", err)
	}

	if gp.config.Verbose {
		fmt.Printf("PDist indices built in %v\n", time.Since(start))
	}

	return nil
}

// GetSupernodeCoordinates generates coordinates for a specific supernode
func (gp *GraphProcessor) GetSupernodeCoordinates(datasetID, supernodeID string, k int) (*CoordinateResult, error) {
	if gp.graph == nil {
		if err := gp.LoadGraph(datasetID); err != nil {
			return nil, err
		}
	}

	if gp.hierarchy == nil {
		if err := gp.LoadHierarchy(datasetID, k); err != nil {
			return nil, err
		}
	}

	// Load indices if needed
	if err := gp.loadIndicesIfNeeded(datasetID); err != nil {
		return nil, fmt.Errorf("failed to load indices: %v", err)
	}

	if gp.config.Verbose {
		fmt.Printf("Generating coordinates for supernode: %s\n", supernodeID)
	}

	start := time.Now()

	// Generate coordinates
	result, err := gp.generateCoordinates(supernodeID)
	if err != nil {
		return nil, fmt.Errorf("failed to generate coordinates: %v", err)
	}

	if gp.config.Verbose {
		fmt.Printf("Coordinates generated in %v\n", time.Since(start))
	}

	return result, nil
}