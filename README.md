# Graph Processor Go Module

A Go implementation of hierarchical graph visualization with Personalized PageRank (PPR) computations, ported from C++ for efficient graph processing and coordinate generation.

## Overview

This module provides three main functionalities:

1. **Build DNPR Indices** (`fpsn` algorithm) - Creates degree-normalized PageRank indices
2. **Build PDist Indices** (`taupush` algorithm) - Creates personalized distance indices  
3. **Generate Coordinates** - On-demand 2D coordinate generation for supernodes using MDS

## Installation

```bash
go mod init your-project
# Copy the graphprocessor package files to your project
```

## Quick Start

```go
package main

import (
    "fmt"
    "log"
    "./graphprocessor" // Adjust import path as needed
)

func main() {
    // Initialize processor
    config := graphprocessor.DefaultConfig()
    config.Verbose = true
    processor := graphprocessor.NewGraphProcessor("/path/to/data", config)

    datasetID := "amazon"
    k := 25

    // Step 1: Build DNPR indices (equivalent to: ./approx_dnppr -f amazon -alg fpsn -k 25 -build 1)
    if err := processor.BuildDNPRIndices(datasetID); err != nil {
        log.Fatalf("Failed to build DNPR indices: %v", err)
    }

    // Step 2: Build PDist indices (equivalent to: ./approx_dnppr -f amazon -alg taupush -k 25 -build 1)
    if err := processor.BuildPDistIndices(datasetID, k); err != nil {
        log.Fatalf("Failed to build PDist indices: %v", err)
    }

    // Step 3: Generate coordinates on-demand (equivalent to: ./approx_dnppr -f amazon -alg taupush -k 25 -ondemand c0_l2_5 -embed 1)
    supernodeID := "c0_l2_5"
    coords, err := processor.GetSupernodeCoordinates(datasetID, supernodeID, k)
    if err != nil {
        log.Fatalf("Failed to generate coordinates: %v", err)
    }

    fmt.Printf("Generated coordinates for %s: %d points\n", supernodeID, len(coords.X))
    fmt.Printf("First point: (%.2f, %.2f), radius: %.2f\n", 
        coords.X[0], coords.Y[0], coords.Radii[0])
}
```

## Directory Structure

Your data directory should be organized as follows:

```
/path/to/data/
├── dataset/
│   ├── amazon.txt              # Graph edges
│   └── amazon_attribute.txt    # Graph attributes (n=nodes, m=edges)
├── louvain/
│   ├── mapping-output/
│   │   └── amazon_25.dat       # Supernode to leaf mapping
│   └── hierachy-output/
│       ├── amazon_25.dat       # Supernode hierarchy
│       └── amazon_25.root      # Root supernodes
├── pr_idx/                     # Generated DNPR indices
└── bwd_idx/                    # Generated backward indices
```

## API Reference

### Core Functions

#### `NewGraphProcessor(basePath string, config *Config) *GraphProcessor`
Creates a new graph processor instance.

#### `BuildDNPRIndices(datasetID string) error`
Builds degree-normalized PageRank indices using the forward push algorithm.
- **Input**: Dataset identifier (e.g., "amazon")
- **Output**: Saves indices to `pr_idx/{datasetID}.dnpr`
- **Equivalent to**: `./approx_dnppr -f {datasetID} -alg fpsn -k {k} -build 1`

#### `BuildPDistIndices(datasetID string, k int) error`
Builds personalized distance indices using the tau-push algorithm.
- **Input**: Dataset identifier and cluster size
- **Output**: Saves indices to `bwd_idx/{datasetID}.*`
- **Equivalent to**: `./approx_dnppr -f {datasetID} -alg taupush -k {k} -build 1`

#### `GetSupernodeCoordinates(datasetID, supernodeID string, k int) (*CoordinateResult, error)`
Generates 2D coordinates for a specific supernode on-demand.
- **Input**: Dataset identifier, supernode ID (e.g., "c0_l2_5"), cluster size
- **Output**: Coordinate result with X, Y positions, radii, and metadata
- **Equivalent to**: `./approx_dnppr -f {datasetID} -alg taupush -k {k} -ondemand {supernodeID} -embed 1`

### Data Structures

#### `Config`
Configuration parameters for the processor:

```go
type Config struct {
    Alpha     float64 // PPR damping factor (default: 0.2)
    Delta     float64 // Precision parameter (auto-calculated)
    PFail     float64 // Failure probability (auto-calculated)
    EpR       float64 // Epsilon parameter (default: 0.5)
    Tau       float64 // Threshold parameter (auto-calculated)
    K         int     // Cluster size (default: 25)
    ThreadNum int     // Number of threads (default: 1)
    Verbose   bool    // Verbose output (default: false)
}
```

#### `CoordinateResult`
Result of coordinate generation:

```go
type CoordinateResult struct {
    X        []float64             // X coordinates
    Y        []float64             // Y coordinates  
    Radii    []float64             // Node radii
    Metadata *SupernodeMetadata    // Additional metadata
}
```

#### `SupernodeMetadata`
Metadata about the generated coordinates:

```go
type SupernodeMetadata struct {
    SupernodeName string     // Supernode identifier
    Level         int        // Hierarchy level
    Children      []string   // Child node/supernode names
    NodeWeights   []int      // Node weights (cluster sizes)
    Degrees       []float64  // Average degrees
    DPRValues     []float64  // Degree PageRank values
    LeafNodes     []int      // Leaf node IDs (level 1 only)
    PPRMatrix     [][]float64 // PPR similarity matrix
    PDistMatrix   [][]float64 // Distance matrix
}
```

### Utility Functions

#### `GetGraphStats() map[string]interface{}`
Returns statistics about the loaded graph and indices.

#### `ValidateDataset(datasetID string, k int) []string`
Validates that all required files exist for a dataset. Returns list of missing files.

#### `SaveCoordinates(result *CoordinateResult, filename string) error`
Saves coordinate results to a JSON file.

#### `LoadCoordinates(filename string) (*CoordinateResult, error)`
Loads coordinate results from a JSON file.

## File Formats

### Graph File Format (`dataset/{datasetID}.txt`)
Edge list format, one edge per line:
```
0 1
0 2
1 3
...
```

### Attribute File Format (`dataset/{datasetID}_attribute.txt`)
```
n=1000
m=5000
```

### Mapping File Format (`louvain/mapping-output/{datasetID}_{k}.dat`)
Supernode to leaf node mapping:
```
c0_l1_0 3 0 1 2
c0_l1_1 2 3 4
...
```

### Hierarchy File Format (`louvain/hierachy-output/{datasetID}_{k}.dat`)
Supernode hierarchy:
```
c0_l2_0 2 0 1
c0_l2_1 3 2 3 4
...
```

## Error Handling

All functions return appropriate errors. Common error scenarios:

- Missing data files
- Corrupted indices
- Invalid supernode IDs
- Memory allocation failures

Always check return values:

```go
if err := processor.BuildDNPRIndices(datasetID); err != nil {
    log.Printf("Error building indices: %v", err)
    return err
}
```

## Performance Notes

- Building indices is CPU-intensive and may take several minutes for large graphs
- Coordinate generation is typically fast (seconds) once indices are built
- Memory usage scales with graph size; ensure sufficient RAM
- Indices are cached on disk for reuse across sessions

## Thread Safety

The GraphProcessor uses mutex locks for thread safety. Multiple goroutines can safely:
- Call read operations simultaneously
- Generate coordinates for different supernodes

Avoid concurrent index building operations on the same processor instance.

## Example: Complete Workflow

```go
package main

import (
    "encoding/json"
    "fmt"
    "log"
    "os"
    "./graphprocessor"
)

func main() {
    // Setup
    config := graphprocessor.DefaultConfig()
    config.Verbose = true
    processor := graphprocessor.NewGraphProcessor("/data", config)

    datasetID := "amazon"
    k := 25

    // Validate dataset
    missing := processor.ValidateDataset(datasetID, k)
    if len(missing) > 0 {
        log.Fatalf("Missing files: %v", missing)
    }

    // Build indices (equivalent to build commands)
    fmt.Println("Building DNPR indices...")
    if err := processor.BuildDNPRIndices(datasetID); err != nil {
        log.Fatalf("DNPR build failed: %v", err)
    }

    fmt.Println("Building PDist indices...")
    if err := processor.BuildPDistIndices(datasetID, k); err != nil {
        log.Fatalf("PDist build failed: %v", err)
    }

    // Generate coordinates on-demand
    fmt.Println("Generating coordinates...")
    coords, err := processor.GetSupernodeCoordinates(datasetID, "c0_l2_0", k)
    if err != nil {
        log.Fatalf("Coordinate generation failed: %v", err)
    }

    // Save results
    if err := processor.SaveCoordinates(coords, "output/coordinates.json"); err != nil {
        log.Printf("Warning: Could not save coordinates: %v", err)
    }

    // Print stats
    stats := processor.GetGraphStats()
    fmt.Println("Graph Statistics:")
    for key, value := range stats {
        fmt.Printf("  %s: %v\n", key, value)
    }

    fmt.Printf("Generated %d coordinate points\n", len(coords.X))
}
```

This completes the documentation for the Graph Processor Go module. The implementation provides a clean, thread-safe API that replicates the functionality of the original C++ code while being idiomatic Go.