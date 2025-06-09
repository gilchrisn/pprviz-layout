// fileio.go - File I/O operations and utility functions

package graphprocessor

import (
	"encoding/gob"
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"math"
	"math/rand"
)

// loadIndicesIfNeeded loads necessary indices if they haven't been loaded
func (gp *GraphProcessor) loadIndicesIfNeeded(datasetID string) error {
	gp.mu.Lock()
	defer gp.mu.Unlock()

	// Load DNPR if not loaded
	if gp.pprIndex == nil || len(gp.pprIndex.DegreePageRank) == 0 {
		dnprPath := filepath.Join(gp.basePath, "pr_idx", datasetID+".dnpr")
		if dnpr, err := gp.loadDNPR(dnprPath); err == nil {
			if gp.pprIndex == nil {
				gp.pprIndex = &PPRIndex{}
			}
			gp.pprIndex.DegreePageRank = dnpr
		} else if gp.config.Verbose {
			fmt.Printf("Could not load DNPR: %v\n", err)
		}
	}

	// Load backward index if not loaded
	if gp.pprIndex == nil || len(gp.pprIndex.BackwardIdxTarget) == 0 {
		bwdPath := filepath.Join(gp.basePath, "bwd_idx", datasetID)
		if err := gp.loadBackwardIndex(bwdPath); err != nil && gp.config.Verbose {
			fmt.Printf("Could not load backward index: %v\n", err)
		}
	}

	return nil
}

// saveDNPR saves degree-normalized PageRank to file
func (gp *GraphProcessor) saveDNPR(filename string, dnpr []float64) error {
	// Ensure directory exists
	if err := os.MkdirAll(filepath.Dir(filename), 0755); err != nil {
		return err
	}

	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	encoder := gob.NewEncoder(file)
	return encoder.Encode(dnpr)
}

// loadDNPR loads degree-normalized PageRank from file
func (gp *GraphProcessor) loadDNPR(filename string) ([]float64, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var dnpr []float64
	decoder := gob.NewDecoder(file)
	err = decoder.Decode(&dnpr)
	return dnpr, err
}

// saveBackwardIndex saves backward push index to files
func (gp *GraphProcessor) saveBackwardIndex(basePath string) error {
	if gp.pprIndex == nil {
		return fmt.Errorf("no PPR index to save")
	}

	// Ensure directory exists
	if err := os.MkdirAll(filepath.Dir(basePath), 0755); err != nil {
		return err
	}

	// Save targets
	targetFile := basePath + ".target"
	if err := gp.saveIntSlice(targetFile, gp.pprIndex.BackwardIdxTarget); err != nil {
		return fmt.Errorf("failed to save targets: %v", err)
	}

	// Save backward index values
	bwdFile := basePath + ".bwdidx"
	if err := gp.saveFloat64Slice(bwdFile, gp.pprIndex.BackwardIdxInCluster); err != nil {
		return fmt.Errorf("failed to save backward index: %v", err)
	}

	// Save info map
	infoFile := basePath + ".bwdidx.info"
	if err := gp.saveBackwardInfo(infoFile, gp.pprIndex.BackwardIdxInfo); err != nil {
		return fmt.Errorf("failed to save backward info: %v", err)
	}

	return nil
}

// loadBackwardIndex loads backward push index from files
func (gp *GraphProcessor) loadBackwardIndex(basePath string) error {
	if gp.pprIndex == nil {
		gp.pprIndex = &PPRIndex{}
	}

	// Load targets
	targetFile := basePath + ".target"
	targets, err := gp.loadIntSlice(targetFile)
	if err != nil {
		return fmt.Errorf("failed to load targets: %v", err)
	}
	gp.pprIndex.BackwardIdxTarget = targets

	// Load backward index values
	bwdFile := basePath + ".bwdidx"
	bwdIdx, err := gp.loadFloat64Slice(bwdFile)
	if err != nil {
		return fmt.Errorf("failed to load backward index: %v", err)
	}
	gp.pprIndex.BackwardIdxInCluster = bwdIdx

	// Load info map
	infoFile := basePath + ".bwdidx.info"
	info, err := gp.loadBackwardInfo(infoFile)
	if err != nil {
		return fmt.Errorf("failed to load backward info: %v", err)
	}
	gp.pprIndex.BackwardIdxInfo = info

	return nil
}

// saveIntSlice saves an integer slice to file
func (gp *GraphProcessor) saveIntSlice(filename string, data []int) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	encoder := gob.NewEncoder(file)
	return encoder.Encode(data)
}

// loadIntSlice loads an integer slice from file
func (gp *GraphProcessor) loadIntSlice(filename string) ([]int, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var data []int
	decoder := gob.NewDecoder(file)
	err = decoder.Decode(&data)
	return data, err
}

// saveFloat64Slice saves a float64 slice to file
func (gp *GraphProcessor) saveFloat64Slice(filename string, data []float64) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	encoder := gob.NewEncoder(file)
	return encoder.Encode(data)
}

// loadFloat64Slice loads a float64 slice from file
func (gp *GraphProcessor) loadFloat64Slice(filename string) ([]float64, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var data []float64
	decoder := gob.NewDecoder(file)
	err = decoder.Decode(&data)
	return data, err
}

// saveBackwardInfo saves backward index info map
func (gp *GraphProcessor) saveBackwardInfo(filename string, info map[int][2]int) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	encoder := gob.NewEncoder(file)
	return encoder.Encode(info)
}

// loadBackwardInfo loads backward index info map
func (gp *GraphProcessor) loadBackwardInfo(filename string) (map[int][2]int, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var info map[int][2]int
	decoder := gob.NewDecoder(file)
	err = decoder.Decode(&info)
	return info, err
}

// SaveCoordinates saves coordinate results to JSON file
func (gp *GraphProcessor) SaveCoordinates(result *CoordinateResult, filename string) error {
	// Ensure directory exists
	if err := os.MkdirAll(filepath.Dir(filename), 0755); err != nil {
		return err
	}

	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	encoder := json.NewEncoder(file)
	encoder.SetIndent("", "  ")
	return encoder.Encode(result)
}

// LoadCoordinates loads coordinate results from JSON file
func (gp *GraphProcessor) LoadCoordinates(filename string) (*CoordinateResult, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var result CoordinateResult
	decoder := json.NewDecoder(file)
	err = decoder.Decode(&result)
	return &result, err
}

// fileExists checks if a file exists
func fileExists(filename string) bool {
	_, err := os.Stat(filename)
	return !os.IsNotExist(err)
}

// ensureDir ensures a directory exists
func ensureDir(dirPath string) error {
	return os.MkdirAll(dirPath, 0755)
}

// Utility functions for algorithms

// randomWalk performs a single random walk from start node
func (gp *GraphProcessor) randomWalk(start int, rng *rand.Rand) int {
	current := start
	
	if len(gp.graph.Edges[start]) == 0 {
		return start
	}

	for {
		if rng.Float64() < gp.config.Alpha {
			return current
		}
		
		neighbors := gp.graph.Edges[current]
		if len(neighbors) == 0 {
			current = start
		} else {
			current = neighbors[rng.Intn(len(neighbors))]
		}
	}
}

// computePPRWithRandomWalk computes PPR using random walks
func (gp *GraphProcessor) computePPRWithRandomWalk(sources []int, insideNodes map[int]bool, omega float64, rsum float64, rng *rand.Rand) map[int]float64 {
	result := make(map[int]float64)
	
	if omega*rsum < 1 {
		return result
	}

	idealWalkNumber := omega * rsum
	increment := rsum / idealWalkNumber
	
	numWalks := int(math.Ceil(idealWalkNumber))
	
	for i := 0; i < numWalks; i++ {
		source := sources[rng.Intn(len(sources))]
		destination := gp.randomWalk(source, rng)
		
		if insideNodes[destination] {
			result[destination] += increment
		}
	}
	
	return result
}

// Helper functions for coordinate generation

// getSupernodeChildren returns children of a supernode
func (gp *GraphProcessor) getSupernodeChildren(supernode string) ([]string, error) {
	_, level := gp.extractSupernodeInfo(supernode)
	
	if level == 1 {
		// Leaf level - return individual nodes
		if leafNodes, exists := gp.hierarchy.Super2Leaf[supernode]; exists {
			children := make([]string, len(leafNodes))
			for i, node := range leafNodes {
				children[i] = fmt.Sprintf("node_%d", node)
			}
			return children, nil
		}
	} else {
		// Higher level - return child supernodes
		if childIDs, exists := gp.hierarchy.Super2Super[supernode]; exists {
			component, _ := gp.extractSupernodeInfo(supernode)
			children := make([]string, len(childIDs))
			for i, childID := range childIDs {
				children[i] = fmt.Sprintf("c%d_l%d_%d", component, level-1, childID)
			}
			return children, nil
		}
	}
	
	return nil, fmt.Errorf("supernode %s not found", supernode)
}

// getLeafNodes returns all leaf nodes for a supernode
func (gp *GraphProcessor) getLeafNodes(supernode string) ([]int, error) {
	if leafNodes, exists := gp.hierarchy.Super2Leaf[supernode]; exists {
		return leafNodes, nil
	}
	return nil, fmt.Errorf("supernode %s not found", supernode)
}

// GetGraphStats returns statistics about the loaded graph
func (gp *GraphProcessor) GetGraphStats() map[string]interface{} {
	stats := make(map[string]interface{})
	
	if gp.graph != nil {
		stats["nodes"] = gp.graph.N
		stats["edges"] = gp.graph.M
		stats["avg_degree"] = gp.graph.DBar
		stats["max_level"] = gp.graph.MaxLevel
	}
	
	if gp.hierarchy != nil {
		stats["supernodes"] = len(gp.hierarchy.Super2Leaf)
		stats["level1_clusters"] = len(gp.hierarchy.Level1Cluster)
		stats["root_supernodes"] = len(gp.hierarchy.Root)
	}
	
	if gp.pprIndex != nil {
		stats["dnpr_loaded"] = len(gp.pprIndex.DegreePageRank) > 0
		stats["backward_index_loaded"] = len(gp.pprIndex.BackwardIdxTarget) > 0
	}
	
	return stats
}

// ValidateDataset checks if all required files exist for a dataset
func (gp *GraphProcessor) ValidateDataset(datasetID string, k int) []string {
	var missing []string
	
	// Check graph files
	graphPath := filepath.Join(gp.basePath, "dataset", datasetID+".txt")
	if !fileExists(graphPath) {
		missing = append(missing, fmt.Sprintf("Graph file: %s", graphPath))
	}
	
	attrPath := filepath.Join(gp.basePath, "dataset", datasetID+"_attribute.txt")
	if !fileExists(attrPath) {
		missing = append(missing, fmt.Sprintf("Attribute file: %s", attrPath))
	}
	
	// Check hierarchy files
	mapPath := filepath.Join(gp.basePath, "louvain", "mapping-output", fmt.Sprintf("%s_%d.dat", datasetID, k))
	if !fileExists(mapPath) {
		missing = append(missing, fmt.Sprintf("Mapping file: %s", mapPath))
	}
	
	hierPath := filepath.Join(gp.basePath, "louvain", "hierachy-output", fmt.Sprintf("%s_%d.dat", datasetID, k))
	if !fileExists(hierPath) {
		missing = append(missing, fmt.Sprintf("Hierarchy file: %s", hierPath))
	}
	
	rootPath := filepath.Join(gp.basePath, "louvain", "hierachy-output", fmt.Sprintf("%s_%d.root", datasetID, k))
	if !fileExists(rootPath) {
		missing = append(missing, fmt.Sprintf("Root file: %s", rootPath))
	}
	
	return missing
}

// SetVerbose enables or disables verbose output
func (gp *GraphProcessor) SetVerbose(verbose bool) {
	gp.config.Verbose = verbose
}

// GetConfig returns a copy of the current configuration
func (gp *GraphProcessor) GetConfig() *Config {
	config := *gp.config
	return &config
}

// UpdateConfig updates the configuration
func (gp *GraphProcessor) UpdateConfig(config *Config) {
	gp.mu.Lock()
	defer gp.mu.Unlock()
	gp.config = config
}