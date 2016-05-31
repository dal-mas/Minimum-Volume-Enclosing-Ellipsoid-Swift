//
//  MinimumVolumeEnclsingEllipsoidSwift.swift
//  
//
//  Created by Felipe Eulalio on 31/05/16.
//
//

import Foundation

/**
Finds the best fit ellipsoid for the given matrix that respects the tollerance
- parameters:
	- P: Matrix that contains N 3D points that are going to be used
	- tollerance: The error tollerance to stop the algorithm
- returns: A matrix containing the values from A, B and C of the parametric equation and the center of the ellipsoid
*/
func findBestFit(P: matrix, tollerance: Double = 0.01) -> ([Double], [Double])
{
	let (size, dimension) = P.shape
	let Q = asarray(P.T.flat.grid + ones(size).grid).reshape((dimension + 1, size)).T
	
	var err = 1 + tollerance
	var u = 1/size * ones(size)
	
	// Khachiyan Algorithm
	while err > tollerance {
		let V = dot(dot(Q.T, y: diag(u)), y: Q)
		let M = diag(dot(dot(Q, y: V.I), y: Q.T))
		let j = argmax(M)
		let max = M[j]
		let stepSize = (max - dimension - 1) / ((dimension + 1) * (max - 1))
		var newU = (1 - stepSize) * u
		
		newU[j] = newU[j] + stepSize
		
		err = norm(newU - u)
		u = newU
	}
	
	let U = diag(u)
	
	// Center of the ellipse
	let center = P.T.dot(u)
	
	let Pu = dot(u.reshape((1, u.count)), y: P)
	let A = (1/dimension) * inv(dot(dot(P.T, y: U), y: P) - (dot(Pu.T, y: Pu)))
	
	// Uses the svd to get the values from 1/(A^2), 1/(B^2) and 1/(C^2).
	// If needed, you can change this code to return the rotation from the 
	// ellipsoid by also returning the rotation
	let (_, values, rotation) = svd(A)
	
	// Get the values of A, B and C
	let radii = 1/sqrt(values)
	
	return (radii.grid, center.grid)
}

/**
Auxiliary method to the Swix Library, that initiates a *ndarray* from the diagonal from a given *matrix*
- parameter M: *matrix* that the main diagonal will be the *ndarray*
- returns: An *ndarray* that contains the main diagonal from the given *matrix*
*/
func diag(M: matrix) -> ndarray
{
	var grid: [Double] = []
	
	for i in 0..<M.rows {
		for j in 0..<M.columns {
			if i == j { grid.append(M[i, j]) }
		}
	}
	
	var array = ndarray(n: grid.count)
	array.grid = grid
	
	return array
}

