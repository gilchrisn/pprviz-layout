// example_test.go - Example usage and basic testing

package main

import (
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strconv"
	"time"
	"./graphprocessor"
)

func main() {
	// Example usage of the GraphProcessor
	
	// Check if data path is provided
	if len(os.Args) < 2 {
		fmt.Println("Usage: go run example_test.go <data_path> [dataset_id] [k]")
		fmt.Println("Example: go run example_test.go /path/to/data amazon 25")
		os.Exit(1)
	}

	basePath := os.Args[1]
	datasetID := "amazon"
	k := 25

	if len(os.Args) > 2 {
		datasetID = os.Args[2]
	}
	if len(os.Args) > 3 {
		if n, err := strconv.Atoi(os.Args[3]); err == nil {
			k = n
		}
	}

	fmt.Printf("Testing GraphProcessor with dataset: %s, k: %d\n", datasetID, k)

	// Create processor with verbose output
	config := graphprocessor.DefaultConfig()
	config.Verbose = true
	processor := graphprocessor.NewGraphProcessor(basePath, config)

	// Test 1: Validate dataset files
	fmt.Println("\n=== Step 1: Validating Dataset ===")
	missing := processor.ValidateDataset(datasetID, k)
	if len(missing) > 0 {
		fmt.Printf("Warning: Missing files:\n")
		for _, file := range missing {
			fmt.Printf("  - %s\n", file)
		}
		fmt.Println("Some operations may fail without these files.")
	} else {
		fmt.Println("All required files found!")
	}

	// Test 2: Load graph and hierarchy
	fmt.Println("\n=== Step 2: Loading Graph and Hierarchy ===")
	if err := processor.LoadGraph(datasetID); err != nil {
		log.Printf("Warning: Could not load graph: %v", err)
	}

	if err := processor.LoadHierarchy(datasetID, k); err != nil {
		log.Printf("Warning: Could not load hierarchy: %v", err)
	}

	// Print graph statistics
	stats := processor.GetGraphStats()
	fmt.Println("Graph Statistics:")
	for key, value := range stats {
		fmt.Printf("  %s: %v\n", key, value)
	}

	// Test 3: Build DNPR indices (equivalent to fpsn algorithm)
	fmt.Println("\n=== Step 3: Building DNPR Indices (fpsn) ===")
	if err := processor.BuildDNPRIndices(datasetID); err != nil {
		log.Printf("Error building DNPR indices: %v", err)
	} else {
		fmt.Println("DNPR indices built successfully!")
	}

	// Test 4: Build PDist indices (equivalent to taupush algorithm)
	fmt.Println("\n=== Step 4: Building PDist Indices (taupush) ===")
	if err := processor.BuildPDistIndices(datasetID, k); err != nil {
		log.Printf("Error building PDist indices: %v", err)
	} else {
		fmt.Println("PDist indices built successfully!")
	}

	// Test 5: Generate coordinates on-demand
	fmt.Println("\n=== Step 5: Generating Coordinates On-Demand ===")
	
	// Try to generate coordinates for a sample supernode
	testSupernodes := []string{"c0_l2_0", "c0_l1_0", "c0_l3_0"}
	
	for _, supernodeID := range testSupernodes {
		fmt.Printf("Attempting to generate coordinates for: %s\n", supernodeID)
		
		coords, err := processor.GetSupernodeCoordinates(datasetID, supernodeID, k)
		if err != nil {
			log.Printf("  Error: %v", err)
			continue
		}

		fmt.Printf("  Success! Generated %d coordinate points\n", len(coords.X))
		fmt.Printf("  Level: %d, Children: %d\n", coords.Metadata.Level, len(coords.Metadata.Children))
		
		if len(coords.X) > 0 {
			fmt.Printf("  First point: (%.2f, %.2f), radius: %.2f\n", 
				coords.X[0], coords.Y[0], coords.Radii[0])
		}

		// Save coordinates to file
		outputFile := fmt.Sprintf("coordinates_%s.json", supernodeID)
		if err := processor.SaveCoordinates(coords, outputFile); err != nil {
			log.Printf("  Warning: Could not save coordinates: %v", err)
		} else {
			fmt.Printf("  Coordinates saved to: %s\n", outputFile)
		}
		
		break // Only test the first successful one
	}

	fmt.Println("\n=== Test Complete ===")
	fmt.Println("GraphProcessor test completed successfully!")
}

// Helper function to create a simple test graph if needed
func createTestGraph(basePath, datasetID string) error {
	// Ensure directories exist
	datasetDir := filepath.Join(basePath, "dataset")
	if err := os.MkdirAll(datasetDir, 0755); err != nil {
		return err
	}

	louvainMappingDir := filepath.Join(basePath, "louvain", "mapping-output")
	if err := os.MkdirAll(louvainMappingDir, 0755); err != nil {
		return err
	}

	louvainHierarchyDir := filepath.Join(basePath, "louvain", "hierachy-output")
	if err := os.MkdirAll(louvainHierarchyDir, 0755); err != nil {
		return err
	}

	// Create simple test graph (triangle)
	graphFile := filepath.Join(datasetDir, datasetID+".txt")
	if err := os.WriteFile(graphFile, []byte("0 1\n1 2\n2 0\n"), 0644); err != nil {
		return err
	}

	// Create attributes file
	attrFile := filepath.Join(datasetDir, datasetID+"_attribute.txt")
	if err := os.WriteFile(attrFile, []byte("n=3\nm=3\n"), 0644); err != nil {
		return err
	}

	// Create simple mapping file
	mapFile := filepath.Join(louvainMappingDir, fmt.Sprintf("%s_25.dat", datasetID))
	mapContent := "c0_l1_0 3 0 1 2\n"
	if err := os.WriteFile(mapFile, []byte(mapContent), 0644); err != nil {
		return err
	}

	// Create simple hierarchy file
	hierFile := filepath.Join(louvainHierarchyDir, fmt.Sprintf("%s_25.dat", datasetID))
	hierContent := "c0_l2_0 1 0\n"
	if err := os.WriteFile(hierFile, []byte(hierContent), 0644); err != nil {
		return err
	}

	// Create root file
	rootFile := filepath.Join(louvainHierarchyDir, fmt.Sprintf("%s_25.root", datasetID))
	if err := os.WriteFile(rootFile, []byte("c0_l2_0\n"), 0644); err != nil {
		return err
	}

	return nil
}

// benchmarkCoordinateGeneration runs a simple benchmark
func benchmarkCoordinateGeneration(processor *graphprocessor.GraphProcessor, datasetID, supernodeID string, k int) {
	fmt.Println("\n=== Benchmark: Coordinate Generation ===")
	
	iterations := 5
	var totalTime time.Duration
	
	for i := 0; i < iterations; i++ {
		start := time.Now()
		_, err := processor.GetSupernodeCoordinates(datasetID, supernodeID, k)
		elapsed := time.Since(start)
		
		if err != nil {
			fmt.Printf("Iteration %d failed: %v\n", i+1, err)
			continue
		}
		
		totalTime += elapsed
		fmt.Printf("Iteration %d: %v\n", i+1, elapsed)
	}
	
	avgTime := totalTime / time.Duration(iterations)
	fmt.Printf("Average time over %d iterations: %v\n", iterations, avgTime)
}