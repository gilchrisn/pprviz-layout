// algorithms.go - Core PPR and coordinate generation algorithms

package graphprocessor

import (
	"fmt"
	"math"
	// "math/rand"
	// "sort"
)

// buildDNPR builds degree-normalized PageRank using power iteration
func (gp *GraphProcessor) buildDNPR(iterations int) ([]float64, error) {
	n := gp.graph.N
	pr := make([]float64, n)
	residuals := make([]float64, n)

	// Initialize residuals with degree-based values
	for i := 0; i < n; i++ {
		residuals[i] = float64(gp.graph.Degrees[i]) / (2.0 * float64(gp.graph.M))
	}

	newResiduals := make([]float64, n)
	rSum := 1.0

	for iter := 0; iter < iterations && rSum > 1e-9; iter++ {
		for id := 0; id < n; id++ {
			degree := gp.graph.Degrees[id]
			if degree == 0 {
				continue
			}

			alphaResidual := gp.config.Alpha * residuals[id]
			pr[id] += alphaResidual
			rSum -= alphaResidual
			increment := (residuals[id] - alphaResidual) / float64(degree)
			residuals[id] = 0

			for _, nid := range gp.graph.Edges[id] {
				newResiduals[nid] += increment
			}
		}
		residuals, newResiduals = newResiduals, residuals
	}

	return pr, nil
}

// buildBackwardPushIndex builds backward push index for high-degree nodes
func (gp *GraphProcessor) buildBackwardPushIndex() error {
	if len(gp.pprIndex.DegreePageRank) == 0 {
		return fmt.Errorf("DNPR not loaded")
	}

	// Create leaf to cluster mapping
	leaf2Cluster := make([]int, gp.graph.N)
	for i, cluster := range gp.hierarchy.Level1Cluster {
		for _, leaf := range cluster {
			if leaf < len(leaf2Cluster) {
				leaf2Cluster[leaf] = i
			}
		}
	}

	tau := gp.config.Tau
	var targets []int
	var infoOffsets []int

	// Collect high-degree targets
	totalSize := 0
	for i := 0; i < gp.graph.N; i++ {
		if gp.pprIndex.DegreePageRank[i] > tau {
			targets = append(targets, i)
			clusterSize := len(gp.hierarchy.Level1Cluster[leaf2Cluster[i]])
			infoOffsets = append(infoOffsets, clusterSize)
			totalSize += clusterSize
		}
	}

	if len(targets) == 0 {
		if gp.config.Verbose {
			fmt.Println("No targets found for backward push")
		}
		return nil
	}

	// Build backward indices
	bwdIdxInCluster := make([]float64, totalSize)
	bwdIdxInfo := make(map[int][2]int)

	offset := 0
	for _, target := range targets {
		cluster := leaf2Cluster[target]
		clusterNodes := gp.hierarchy.Level1Cluster[cluster]
		
		// Compute backward push for this target
		bwdResult, err := gp.backwardPush(target, tau)
		if err != nil {
			continue
		}

		// Store results for cluster nodes
		for j, node := range clusterNodes {
			if reserve, exists := bwdResult.Reserve[node]; exists {
				bwdIdxInCluster[offset+j] = float64(gp.graph.Degrees[node]) * reserve
			}
		}

		bwdIdxInfo[target] = [2]int{offset, len(clusterNodes)}
		offset += len(clusterNodes)
	}

	// Store in PPR index
	if gp.pprIndex == nil {
		gp.pprIndex = &PPRIndex{}
	}
	gp.pprIndex.BackwardIdxTarget = targets
	gp.pprIndex.BackwardIdxInCluster = bwdIdxInCluster
	gp.pprIndex.BackwardIdxInfo = bwdIdxInfo

	if gp.config.Verbose {
		fmt.Printf("Built backward push index for %d targets\n", len(targets))
	}

	return nil
}

// backwardPush performs backward push from a target node
func (gp *GraphProcessor) backwardPush(target int, rmax float64) (*ForwardIndex, error) {
	reserve := make(map[int]float64)
	residual := make(map[int]float64)
	
	queue := []int{target}
	residual[target] = 1.0

	for len(queue) > 0 {
		v := queue[0]
		queue = queue[1:]

		vResidue, exists := residual[v]
		if !exists || vResidue == 0 {
			continue
		}

		residual[v] = 0
		reserve[v] += vResidue * gp.config.Alpha

		remaining := (1.0 - gp.config.Alpha) * vResidue

		for _, next := range gp.graph.Edges[v] {
			if gp.graph.Degrees[next] == 0 {
				continue
			}

			oldResidual := residual[next]
			newResidual := remaining / float64(gp.graph.Degrees[next])
			residual[next] = oldResidual + newResidual

			if oldResidual/float64(gp.graph.Degrees[next]) <= rmax && 
			   residual[next]/float64(gp.graph.Degrees[next]) > rmax {
				queue = append(queue, next)
			}
		}
	}

	return &ForwardIndex{Reserve: reserve, Residual: residual}, nil
}

// forwardPush performs forward push from source nodes
func (gp *GraphProcessor) forwardPush(sources []int, rmax float64) (*ForwardIndex, error) {
	reserve := make(map[int]float64)
	residual := make(map[int]float64)
	
	// Initialize with sources
	for _, source := range sources {
		residual[source] = float64(gp.graph.Degrees[source])
	}

	queue := make([]int, len(sources))
	copy(queue, sources)

	for len(queue) > 0 {
		v := queue[0]
		queue = queue[1:]

		vResidue, exists := residual[v]
		if !exists || vResidue == 0 || gp.graph.Degrees[v] == 0 {
			continue
		}

		if vResidue/float64(gp.graph.Degrees[v]) < rmax {
			continue
		}

		residual[v] = 0
		reserve[v] += vResidue * gp.config.Alpha

		outNeighbors := len(gp.graph.Edges[v])
		if outNeighbors == 0 {
			residual[sources[0]] += vResidue * (1 - gp.config.Alpha)
			continue
		}

		avgPushResidual := ((1.0 - gp.config.Alpha) * vResidue) / float64(outNeighbors)

		for _, next := range gp.graph.Edges[v] {
			oldResidual := residual[next]
			residual[next] = oldResidual + avgPushResidual

			if gp.graph.Degrees[next] > 0 &&
			   oldResidual/float64(gp.graph.Degrees[next]) <= rmax &&
			   residual[next]/float64(gp.graph.Degrees[next]) > rmax {
				queue = append(queue, next)
			}
		}
	}

	return &ForwardIndex{Reserve: reserve, Residual: residual}, nil
}

// computeSuperPPR computes PPR between supernodes using tau-push algorithm
func (gp *GraphProcessor) computeSuperPPR(supernode string) ([][]float64, []string, error) {
	_, level := gp.extractSupernodeInfo(supernode)
	
	var children []string
	var nodeWeights []int

	if level == 1 {
		// Leaf level - use individual nodes
		leafNodes := gp.hierarchy.Super2Leaf[supernode]
		for _, leaf := range leafNodes {
			children = append(children, fmt.Sprintf("node_%d", leaf))
		}
		nodeWeights = make([]int, len(leafNodes))
		for i := range nodeWeights {
			nodeWeights[i] = 1
		}
	} else {
		// Higher level - use child supernodes
		childIDs := gp.hierarchy.Super2Super[supernode]
		component, _ := gp.extractSupernodeInfo(supernode)
		
		for _, childID := range childIDs {
			childName := fmt.Sprintf("c%d_l%d_%d", component, level-1, childID)
			children = append(children, childName)
			if leafNodes, exists := gp.hierarchy.Super2Leaf[childName]; exists {
				nodeWeights = append(nodeWeights, len(leafNodes))
			} else {
				nodeWeights = append(nodeWeights, 1)
			}
		}
	}

	n := len(children)
	if n == 0 {
		return nil, nil, fmt.Errorf("no children found for supernode %s", supernode)
	}

	pprMatrix := make([][]float64, n)
	for i := range pprMatrix {
		pprMatrix[i] = make([]float64, n)
	}

	// Compute PPR for each source
	for i := 0; i < n; i++ {
		var sources []int
		
		if level == 1 {
			// For leaf level, use individual nodes
			leafNodes := gp.hierarchy.Super2Leaf[supernode]
			if i < len(leafNodes) {
				sources = []int{leafNodes[i]}
			}
		} else {
			// For higher levels, use all nodes in the child supernode
			childName := children[i]
			if leafNodes, exists := gp.hierarchy.Super2Leaf[childName]; exists {
				sources = leafNodes
			}
		}

		if len(sources) == 0 {
			continue
		}

		// Perform forward push
		rmax := gp.computeRMax(len(sources))
		fwdResult, err := gp.forwardPush(sources, rmax)
		if err != nil {
			continue
		}

		// Aggregate results for target supernodes
		for j := 0; j < n; j++ {
			var targetNodes []int
			
			if level == 1 {
				leafNodes := gp.hierarchy.Super2Leaf[supernode]
				if j < len(leafNodes) {
					targetNodes = []int{leafNodes[j]}
				}
			} else {
				targetName := children[j]
				if leafNodes, exists := gp.hierarchy.Super2Leaf[targetName]; exists {
					targetNodes = leafNodes
				}
			}

			totalReserve := 0.0
			for _, target := range targetNodes {
				if reserve, exists := fwdResult.Reserve[target]; exists {
					totalReserve += reserve
				}
			}

			if len(targetNodes) > 0 {
				pprMatrix[i][j] = totalReserve / float64(len(targetNodes))
			}
		}
	}

	return pprMatrix, children, nil
}

// computeRMax computes the rmax parameter for forward push
func (gp *GraphProcessor) computeRMax(sourceSize int) float64 {
	deltap := gp.graph.DBar * gp.config.Delta
	maxVal := 1.0 - math.Pow(float64(gp.graph.N), -2*gp.config.EpR)
	minVal := 1.0 - math.Exp(-2*gp.config.EpR)
	
	epsilonP := math.Min(maxVal, math.Max(1.0-math.Pow(2*deltap/math.E, gp.config.EpR), minVal))
	
	return epsilonP * deltap / gp.config.Tau / float64(gp.graph.M) / 2.0
}

// generateCoordinates generates 2D coordinates using MDS
func (gp *GraphProcessor) generateCoordinates(supernodeID string) (*CoordinateResult, error) {
	// Compute PPR matrix
	pprMatrix, children, err := gp.computeSuperPPR(supernodeID)
	if err != nil {
		return nil, err
	}

	if len(children) == 0 {
		return nil, fmt.Errorf("no children found for supernode %s", supernodeID)
	}

	n := len(children)
	
	// Convert PPR to distance matrix
	distMatrix := gp.pprToDistance(pprMatrix)
	
	// Compute node weights
	nodeWeights := make([]int, n)
	degrees := make([]float64, n)
	dprValues := make([]float64, n)
	
	_, level := gp.extractSupernodeInfo(supernodeID)
	var leafNodes []int

	for i, child := range children {
		if level == 1 {
			// Individual leaf nodes
			leafNodes := gp.hierarchy.Super2Leaf[supernodeID]
			if i < len(leafNodes) {
				leaf := leafNodes[i]
				nodeWeights[i] = 1
				degrees[i] = float64(gp.graph.Degrees[leaf])
				if i < len(gp.pprIndex.DegreePageRank) {
					dprValues[i] = gp.pprIndex.DegreePageRank[leaf]
				}
				leafNodes = append(leafNodes, leaf)
			}
		} else {
			// Child supernodes
			if childLeaves, exists := gp.hierarchy.Super2Leaf[child]; exists {
				nodeWeights[i] = len(childLeaves)
				
				totalDegree := 0.0
				totalDPR := 0.0
				for _, leaf := range childLeaves {
					totalDegree += float64(gp.graph.Degrees[leaf])
					if leaf < len(gp.pprIndex.DegreePageRank) {
						totalDPR += gp.pprIndex.DegreePageRank[leaf]
					}
				}
				degrees[i] = totalDegree / float64(len(childLeaves))
				dprValues[i] = totalDPR / float64(len(childLeaves))
			} else {
				nodeWeights[i] = 1
				degrees[i] = 1.0
				dprValues[i] = 0.0
			}
		}
	}

	// Add radii to distance matrix
	radii := gp.computeRadii(nodeWeights, distMatrix)
	distMatrixWithRadii := gp.addRadiiToDistances(distMatrix, radii)

	// Perform MDS
	coordinates, err := gp.performMDS(distMatrixWithRadii)
	if err != nil {
		return nil, err
	}

	// Create metadata
	metadata := &SupernodeMetadata{
		SupernodeName: supernodeID,
		Level:         level,
		Children:      children,
		NodeWeights:   nodeWeights,
		Degrees:       degrees,
		DPRValues:     dprValues,
		PPRMatrix:     pprMatrix,
		PDistMatrix:   distMatrix,
	}

	if level == 1 && len(leafNodes) > 0 {
		metadata.LeafNodes = leafNodes
	}

	return &CoordinateResult{
		X:        coordinates[0],
		Y:        coordinates[1],
		Radii:    radii,
		Metadata: metadata,
	}, nil
}

// pprToDistance converts PPR values to distance values
func (gp *GraphProcessor) pprToDistance(pprMatrix [][]float64) [][]float64 {
	n := len(pprMatrix)
	distMatrix := make([][]float64, n)
	for i := range distMatrix {
		distMatrix[i] = make([]float64, n)
	}

	maxVal := 2 * math.Log2(float64(gp.graph.N))
	minVal := 2.0

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			val := pprMatrix[i][j] + pprMatrix[j][i]
			if val == 0 {
				val = maxVal
			} else {
				val = 1 - math.Log2(val)
			}
			val = math.Min(math.Max(minVal, val), maxVal)
			distMatrix[i][j] = val
		}
	}

	return distMatrix
}

// computeRadii computes radii for nodes based on their weights
func (gp *GraphProcessor) computeRadii(nodeWeights []int, distMatrix [][]float64) []float64 {
	n := len(nodeWeights)
	radii := make([]float64, n)

	avgSize := 0.0
	avgDist := 0.0
	
	for _, weight := range nodeWeights {
		avgSize += float64(weight)
	}
	avgSize /= float64(n)

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			avgDist += distMatrix[i][j]
		}
	}
	avgDist /= float64(n * n)

	scale := 0.03 * avgDist / math.Sqrt(avgSize)
	
	for i := 0; i < n; i++ {
		radii[i] = scale * math.Sqrt(float64(nodeWeights[i]))
	}

	return radii
}

// addRadiiToDistances adds radii to the distance matrix
func (gp *GraphProcessor) addRadiiToDistances(distMatrix [][]float64, radii []float64) [][]float64 {
	n := len(distMatrix)
	result := make([][]float64, n)
	for i := range result {
		result[i] = make([]float64, n)
		for j := range result[i] {
			result[i][j] = distMatrix[i][j] + radii[i] + radii[j]
		}
	}
	return result
}

// performMDS performs Multidimensional Scaling to generate 2D coordinates
func (gp *GraphProcessor) performMDS(distMatrix [][]float64) ([][]float64, error) {
	n := len(distMatrix)
	if n == 0 {
		return nil, fmt.Errorf("empty distance matrix")
	}

	// Initialize coordinates randomly in a circle
	coordinates := make([][]float64, 2)
	coordinates[0] = make([]float64, n) // X coordinates
	coordinates[1] = make([]float64, n) // Y coordinates

	for i := 0; i < n; i++ {
		angle := 2.0 * math.Pi * float64(i) / float64(n)
		coordinates[0][i] = 50.0 * math.Cos(angle)
		coordinates[1][i] = 50.0 * math.Sin(angle)
	}

	// SMACOF algorithm implementation
	weights := gp.computeWeights(distMatrix)
	
	maxIter := 300
	tolerance := 1e-3
	
	for iter := 0; iter < maxIter; iter++ {
		// Compute current distances
		currentDists := gp.computeCurrentDistances(coordinates)
		
		// Compute stress
		stress := gp.computeStress(weights, distMatrix, currentDists)
		
		// Update coordinates
		newCoordinates := gp.updateCoordinatesSMACOF(coordinates, weights, distMatrix, currentDists)
		
		// Compute new stress
		newCurrentDists := gp.computeCurrentDistances(newCoordinates)
		newStress := gp.computeStress(weights, distMatrix, newCurrentDists)
		
		// Check convergence
		if math.Abs(stress-newStress) < tolerance {
			break
		}
		
		coordinates = newCoordinates
	}

	return coordinates, nil
}

// computeWeights computes weights for MDS
func (gp *GraphProcessor) computeWeights(distMatrix [][]float64) [][]float64 {
	n := len(distMatrix)
	weights := make([][]float64, n)
	for i := range weights {
		weights[i] = make([]float64, n)
		for j := range weights[i] {
			if distMatrix[i][j] > 0 {
				weights[i][j] = 1.0 / (distMatrix[i][j] * distMatrix[i][j])
			}
		}
	}
	return weights
}

// computeCurrentDistances computes current Euclidean distances
func (gp *GraphProcessor) computeCurrentDistances(coordinates [][]float64) [][]float64 {
	n := len(coordinates[0])
	dists := make([][]float64, n)
	for i := range dists {
		dists[i] = make([]float64, n)
		for j := range dists[i] {
			dx := coordinates[0][i] - coordinates[0][j]
			dy := coordinates[1][i] - coordinates[1][j]
			dists[i][j] = math.Sqrt(dx*dx + dy*dy)
		}
	}
	return dists
}

// computeStress computes the stress function for MDS
func (gp *GraphProcessor) computeStress(weights, targetDists, currentDists [][]float64) float64 {
	stress := 0.0
	totalWeight := 0.0
	
	n := len(weights)
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			diff := targetDists[i][j] - currentDists[i][j]
			stress += weights[i][j] * diff * diff
			totalWeight += weights[i][j]
		}
	}
	
	if totalWeight > 0 {
		stress /= totalWeight
	}
	
	return stress
}

// updateCoordinatesSMACOF updates coordinates using SMACOF algorithm
func (gp *GraphProcessor) updateCoordinatesSMACOF(coordinates, weights, targetDists, currentDists [][]float64) [][]float64 {
	n := len(coordinates[0])
	newCoords := make([][]float64, 2)
	newCoords[0] = make([]float64, n)
	newCoords[1] = make([]float64, n)

	// Compute B matrix and update coordinates
	for dim := 0; dim < 2; dim++ {
		for i := 0; i < n; i++ {
			numerator := 0.0
			denominator := 0.0
			
			for j := 0; j < n; j++ {
				if i == j {
					continue
				}
				
				wij := weights[i][j]
				if currentDists[i][j] > 0 {
					bij := wij * targetDists[i][j] / currentDists[i][j]
					numerator += bij * coordinates[dim][j]
					denominator += wij
				}
			}
			
			if denominator > 0 {
				newCoords[dim][i] = numerator / denominator
			} else {
				newCoords[dim][i] = coordinates[dim][i]
			}
		}
	}

	return newCoords
}