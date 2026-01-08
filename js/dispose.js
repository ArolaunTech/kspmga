// Script to dispose of THREE.js node
// Thanks to https://stackoverflow.com/questions/33152132/three-js-collada-whats-the-proper-way-to-dispose-and-release-memory-garbag/33199591#33199591
import * as THREE from 'three';

export function disposeMaterial(material) {
	if (material.map) material.map.dispose();
	if (material.lightMap) material.lightMap.dispose();
	if (material.bumpMap) material.bumpMap.dispose();
	if (material.normalMap) material.normalMap.dispose();
	if (material.specularMap) material.specularMap.dispose();
	if (material.envMap) material.envMap.dispose();
	if (material.alphaMap) material.alphaMap.dispose();
	if (material.aoMap) material.aoMap.dispose();
	if (material.displacementMap) material.displacementMap.dispose();
	if (material.emissiveMap) material.emissiveMap.dispose();
	if (material.gradientMap) material.gradientMap.dispose();
	if (material.metalnessMap) material.metalnessMap.dispose();
	if (material.roughnessMap) material.roughnessMap.dispose();

	material.dispose();
}

export function disposeNode(node) {
	if (node.geometry) node.geometry.dispose();

	if (node.material) {
		if (Array.isArray(node.material)) {
			node.material.forEach(disposeMaterial);
		} else {
			disposeMaterial(node.material);
		}
	}
}

export function disposeHierarchy(node) {
	for (let i = node.children.length - 1; i >= 0; i--) {
		let child = node.children[i];
		disposeHierarchy(child);
		disposeNode(child);
	}
}