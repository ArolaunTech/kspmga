class Vector3 {
	constructor(x, y, z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}

	get norm() {
		return Math.hypot(x, y, z);
	}
};

function dot3(a, b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

function cross(a, b) {
	return new Vector3(
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x
	);
}

function normalize(a) {
	let norm = a.norm;

	return new Vector3(a.x / norm, a.y / norm, a.z / norm);
}