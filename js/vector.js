export class Vector3 {
	constructor(x, y, z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}

	get norm() {
		return Math.hypot(this.x, this.y, this.z);
	}
};

// TO-DO: turn these functions into static functions

export function add3(a, b) {
	return new Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
}

export function dot3(a, b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

export function cross(a, b) {
	return new Vector3(
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x
	);
}

export function mult3(a, b) {
	return new Vector3(a.x * b, a.y * b, a.z * b);
}

export function normalize3(a) {
	return mult3(a, 1 / a.norm);
}