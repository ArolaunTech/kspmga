export function colorArrToHex(col) {
	let r = Math.floor(255 * col[0]);
	let g = Math.floor(255 * col[1]);
	let b = Math.floor(255 * col[2]);

	return 65536 * r + 256 * g + b;
}