export class Complex {
	constructor(a, b) {
		this.a = a;
		this.b = b;
	}

	get sqnorm() {
		return this.a * this.a + this.b * this.b;
	}

	get absvalue() {
		return Math.sqrt(this.sqnorm);
	}

	static fromReal(x) {
		return new Complex(x, 0);
	}

	static add(x, y) {
		return new Complex(x.a + y.a, x.b + y.b);
	}

	static sub(x, y) {
		return new Complex(x.a - y.a, x.b - y.b);
	}

	static mult(x, y) {
		return new Complex(x.a * y.a - x.b * y.b, x.b * y.a + x.a * y.b);
	}

	static multReal(x, y) {
		return new Complex(x.a * y, x.b * y);
	}

	static conjugate(x) {
		return new Complex(x.a, -x.b);
	}

	static reciprocal(x) {
		return new Complex(x.a / x.sqnorm, -x.b / x.sqnorm);
	}

	static divide(x, y) {
		return Complex.mult(x, Complex.reciprocal(y));
	}

	static powInt(x, y) {
		let out = new Complex(1, 0);

		for (let i = 0; i < y; i++) {
			out = Complex.mult(out, x);
		}

		return out;
	}

	static powReal(x, y) {
		if (x.a === 0 && x.b === 0) {
			return x;
		}

		let angle = Math.atan2(x.b, x.a) * y;
		let norm = Math.pow(x.absvalue, y);

		return new Complex(
			norm * Math.cos(angle),
			norm * Math.sin(angle)
		);
	}

	static sqrt(x) {
		if (x.a === 0 && x.b === 0) {
			return x;
		}
		if (x.b === 0 && x.a < 0) {
			return new Complex(0, Math.sqrt(-x.a));
		}
		if (x.b === 0) {
			return new Complex(Math.sqrt(x.a), 0);
		}

		let halfa = (x.a + x.absvalue) / 2;
		let halfb = x.b / 2;

		let sqrtnorm = Math.sqrt(x.absvalue);
		let halfnorm = Math.hypot(halfa, halfb);

		return new Complex(halfa * sqrtnorm / halfnorm, halfb * sqrtnorm / halfnorm);
	}

	static cbrt(x) {
		return Complex.powReal(x, 1/3);
	}
};

export function solveQuadratic(b, c) {
	// Solves the equation
	//  2
	// x  + bx + c = 0
	//
	// analytically.

	let discriminant = b * b - 4 * c;

	if (discriminant <= 0) {
		return [
			new Complex(-b / 2, Math.sqrt(-discriminant) / 2),
			new Complex(-b / 2, -Math.sqrt(-discriminant) / 2)
		];
	}

	let root1;
	if (b < 0) {
		root1 = (-b + Math.sqrt(discriminant)) / 2;
	} else {
		root1 = (-b - Math.sqrt(discriminant)) / 2;
	}

	return [
		new Complex(root1, 0),
		new Complex(c / root1, 0)
	];
}

export function solveQuarticBairstow(p, q, r, s) {
	// Solves the equation
	//  4     3     2
	// x  + px  + qx  + rx + s = 0

	let converged = false;
	let iters = 0;

	let u, v;

	while (!converged && iters < 1000) {
		iters++;
		u = Math.random();
		v = Math.random();

		for (let i = 0; i < 50; i++) {
			let b1 = p - u;
			let b0 = q - v - u * b1;

			let c = r - v * b1 - u * b0;
			let d = s - v * b0;

			let g = b1 - u;
			let h = b0 - v;

			let quotient = 1 / (v * g * g + h * (h - u * g));
			let du = -quotient * (-h * c + g * d);
			let dv = -quotient * (-g * v * c + (g * u - h) * d);

			u += du;
			v += dv;

			if ((c * c + d * d) < 1e-30) {
				converged = true;
				break;
			}
		}
	}

	let b1 = p - u;
	let b0 = q - v - u * b1;

	return solveQuadratic(u, v).concat(solveQuadratic(b1, b0));
}

export function solveQuarticMine(q3, q2, q1, q0) {
	// I could do a scaling step but most of the time it won't be necessary

	// Find coefficients of associated cubic
	let c0 = 4 * q0 * q2 - q0 * q3 * q3 - q1 * q1;
	let c1 = q1 * q3 - 4 * q0;
	let c2 = -q2;

	// Find coefficients of depressed form
	let d0 = 8 * q0 * q2 / 3 - q0 * q3 * q3 - q1 * q1 - 2 * q2 * q2 * q2 / 27 + q1 * q2 * q3 / 3;
	let d1 = q1 * q3 - 4 * q0 - q2 * q2 / 3;

	// Find initial guess for Newton's method
	let xg = 2 * Math.max(Math.abs(c2), Math.sqrt(Math.abs(c1)), Math.cbrt(Math.abs(c0)));
	if (d1 > 0 && d0 > 0) xg = -xg;
	if (d1 < 0 && d0 > 0) {
		let xmin = (-c2 + Math.sqrt(c2 * c2 - 3 * c1)) / 3;
		let ymin = ((xmin + c2) * xmin + c1) * xmin + c0;

		if (ymin > 0) xg = -xg;
	}

	// Solve cubic using Newton's method
	for (let i = 0; i < 50; i++) {
		let y = ((xg + c2) * xg + c1) * xg + c0;
		let dy = (3 * xg + 2 * c2) * xg + c1;

		let dx = y / dy;

		let xgprev = xg;

		xg -= dx;

		if (Math.abs(xg - xgprev) < 1e-12) break;
	}

	// Find u, v from xg
	let discriminant = xg * xg - 4 * q0;
	if (discriminant < 0) discriminant = 0;

	let v = (xg >= 0) ? (xg + Math.sqrt(discriminant)) : (xg - Math.sqrt(discriminant));
	v /= 2;

	let u = v * (q3 * v - q1) / (v * v - q0);

	// Bairstow iterations to refine u, v
	for (let i = 0; i < 15; i++) {
		let b1 = q3 - u;
		let b0 = q2 - v - u * b1;

		let c = q1 - v * b1 - u * b0;
		let d = q0 - v * b0;

		let g = b1 - u;
		let h = b0 - v;

		let quotient = 1 / (v * g * g + h * (h - u * g));
		let du = -quotient * (-h * c + g * d);
		let dv = -quotient * (-g * v * c + (g * u - h) * d);

		let pu = u;
		let pv = v;

		u += du;
		v += dv;

		if (u == pu && v == pv) break;
	}

	// Find other quadratic factor
	let b1 = q3 - u;
	let b0 = q2 - v - u * b1;

	// Our quartic factors into
	//  4     3     2              / 2        \  / 2           \
	// x + q x + q x + q x + q  = | x + ux + v || x + b x + b   |
	//      3     2     1     0    \          /  \     1     0 /

	// Once again no scaling step
	return solveQuadratic(u, v).concat(solveQuadratic(b1, b0));
}

export function realRootsQuartic(p, q, r, s) {
	let roots = solveQuarticMine(p, q, r, s);
	let out = [];

	for (let i = 0; i < 4; i++) {
		//console.log(roots[i]);

		if (Math.abs(roots[i].b) < 1e-10) {
			out.push(roots[i].a);
		}
	}

	out.sort();

	return out;
}