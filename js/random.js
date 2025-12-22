export function normal() {
	// Simple Box-Muller transform
	return Math.sqrt(-2 * Math.log(Math.random())) * Math.cos(2 * Math.PI * Math.random());
}