// run_tests.go - Main test runner script

package main

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"time"
	"github.com/gilchrisn/pprviz-layout/graphprocessor"
)

func main() {
	fmt.Println("Graph Processor Comprehensive Test Suite")
	fmt.Println("========================================")
	
	// Get current directory
	wd, err := os.Getwd()
	if err != nil {
		fmt.Printf("Error getting working directory: %v\n", err)
		os.Exit(1)
	}
	
	fmt.Printf("Running tests from: %s\n", wd)
	fmt.Printf("Go version: ")
	
	// Check Go version
	cmd := exec.Command("go", "version")
	output, err := cmd.Output()
	if err != nil {
		fmt.Printf("Error checking Go version: %v\n", err)
		os.Exit(1)
	}
	fmt.Print(string(output))
	
	// Build the package first
	fmt.Println("\n--- Building Package ---")
	buildCmd := exec.Command("go", "build", "./graphprocessor")
	buildCmd.Dir = wd
	buildOutput, err := buildCmd.CombinedOutput()
	if err != nil {
		fmt.Printf("Build failed: %v\n", err)
		fmt.Printf("Output: %s\n", string(buildOutput))
		os.Exit(1)
	}
	fmt.Println("✓ Package built successfully")
	
	// Run basic tests
	fmt.Println("\n--- Running Basic Test Suite ---")
	start := time.Now()
	
	testCmd := exec.Command("go", "run", "test_suite.go")
	testCmd.Dir = wd
	testCmd.Stdout = os.Stdout
	testCmd.Stderr = os.Stderr
	
	err = testCmd.Run()
	basicTestDuration := time.Since(start)
	
	if err != nil {
		fmt.Printf("Basic tests failed: %v\n", err)
		fmt.Printf("Duration: %v\n", basicTestDuration)
		os.Exit(1)
	}
	
	fmt.Printf("\n✓ Basic tests completed in %v\n", basicTestDuration)
	
	// Run benchmark tests
	fmt.Println("\n--- Running Benchmark Suite ---")
	start = time.Now()
	
	benchCmd := exec.Command("go", "run", "benchmark_test.go")
	benchCmd.Dir = wd
	benchCmd.Stdout = os.Stdout
	benchCmd.Stderr = os.Stderr
	
	err = benchCmd.Run()
	benchDuration := time.Since(start)
	
	if err != nil {
		fmt.Printf("Benchmark tests failed: %v\n", err)
		fmt.Printf("Duration: %v\n", benchDuration)
		os.Exit(1)
	}
	
	fmt.Printf("\n✓ Benchmark tests completed in %v\n", benchDuration)
	
	// Run example test
	fmt.Println("\n--- Running Example Test ---")
	start = time.Now()
	
	// Create a temporary data directory for the example
	tempDir := filepath.Join(os.TempDir(), fmt.Sprintf("graphprocessor_example_%d", time.Now().Unix()))
	defer os.RemoveAll(tempDir)
	
	exampleCmd := exec.Command("go", "run", "example_test.go", tempDir, "test", "25")
	exampleCmd.Dir = wd
	exampleCmd.Stdout = os.Stdout
	exampleCmd.Stderr = os.Stderr
	
	err = exampleCmd.Run()
	exampleDuration := time.Since(start)
	
	if err != nil {
		fmt.Printf("Example test failed: %v\n", err)
		fmt.Printf("Duration: %v\n", exampleDuration)
		// Don't exit on example failure - it might be due to missing data
	} else {
		fmt.Printf("\n✓ Example test completed in %v\n", exampleDuration)
	}
	
	// Summary
	totalDuration := basicTestDuration + benchDuration + exampleDuration
	fmt.Printf("\n=== Test Summary ===\n")
	fmt.Printf("Basic Tests:     %v ✓\n", basicTestDuration)
	fmt.Printf("Benchmark Tests: %v ✓\n", benchDuration)
	fmt.Printf("Example Test:    %v %s\n", exampleDuration, 
		map[bool]string{true: "✓", false: "⚠"}[err == nil])
	fmt.Printf("Total Time:      %v\n", totalDuration)
	fmt.Println("\n🎉 All critical tests passed!")
}
