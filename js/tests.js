export function runTest(system) {
	const start = performance.now();
	const runs = 10000;

	for (let i = 0; i < runs; i++) {
		console.log(system.getState("Kerbin", 1e8 * Math.random()));
	}

	console.log((performance.now() - start) / runs);
}