export class Vector3 {
	constructor(x, y, z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}

	get norm() {
		return Math.hypot(this.x, this.y, this.z);
	}

	static add(a, b) {
		return new Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
	}

	static sub(a, b) {
		return new Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
	}

	static dot(a, b) {
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}

	static cross(a, b) {
		return new Vector3(
			a.y * b.z - a.z * b.y,
			a.z * b.x - a.x * b.z,
			a.x * b.y - a.y * b.x
		);
	}

	static mult(a, b) {
		return new Vector3(a.x * b, a.y * b, a.z * b);
	}

	static normalize(a) {
		return Vector3.mult(a, 1 / a.norm);
	}
};