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

export function solveCubicAnalytic(b, c, d) {
	// Equation: x^3 + bx^2 + cx + d = 0

	let delta0 = b * b - 3 * c;
	let delta1 = 2 * b * b * b - 9 * b * c + 27 * d;

	if (delta0 === 0 && delta1 === 0) {
		return [
			Complex.fromReal(-b / 3),
			Complex.fromReal(-b / 3),
			Complex.fromReal(-b / 3)
		];
	}

	let insideroot = delta1 * delta1 - 4 * delta0 * delta0 * delta0;

	let C;

	if (insideroot <= 0) {
		C = Complex.cbrt(new Complex(delta1 / 2, 0.5 * Math.sqrt(-insideroot)));
	} else if (delta1 >= 0) {
		C = Complex.cbrt(
			Complex.fromReal(
				0.5 * (delta1 + Math.sqrt(insideroot))
			)
		);
	} else {
		C = Complex.cbrt(
			Complex.fromReal(
				0.5 * (delta1 - Math.sqrt(insideroot))
			)
		);
	}

	let mult1 = new Complex(-0.5, Math.sqrt(3) / 2);
	let mult2 = new Complex(-0.5, -Math.sqrt(3) / 2);

	let C0 = C;
	let C1 = Complex.mult(C, mult1);
	let C2 = Complex.mult(C, mult2);

	return [
		Complex.multReal(
			Complex.add(
				Complex.fromReal(b), 
				Complex.add(
					C0, 
					Complex.multReal(Complex.reciprocal(C0), delta0))
			), 
			-1/3
		),
		Complex.multReal(
			Complex.add(
				Complex.fromReal(b), 
				Complex.add(
					C1, 
					Complex.multReal(Complex.reciprocal(C1), delta0))
				),
			-1/3
		),
		Complex.multReal(
			Complex.add(
				Complex.fromReal(b), 
				Complex.add(
					C2, 
					Complex.multReal(Complex.reciprocal(C2), delta0))
				),
			-1/3
		)
	];
}

export function solveQuarticDurandKerner(p, q, r, s) {
	// Solves the equation
	//  4     3     2
	// x  + px  + qx  + rx + s = 0

	let sols = [new Complex(0.4, 0.9), new Complex(-0.65, 0.72), new Complex(-0.908, -0.297), new Complex(-0.0959, -0.936)];

	let maxiters = 100;
	for (let i = 0; i < maxiters; i++) {
		let newsols = structuredClone(sols);
		let maxf = 0;

		for (let idx = 0; idx < 4; idx++) {
			let f = Complex.add(
				Complex.add(
					Complex.add(
						Complex.powInt(sols[idx], 4), 
						Complex.multReal(Complex.powInt(sols[idx], 3), p)
					),
					Complex.add(
						Complex.multReal(Complex.powInt(sols[idx], 2), q),
						Complex.multReal(sols[idx], r)
					)
				),
				Complex.fromReal(s)
			);

			for (let idx2 = 0; idx2 < 4; idx2++) {
				if (idx == idx2) continue;

				f = Complex.divide(f, Complex.sub(sols[idx], sols[idx2]));
			}

			if (f.sqnorm > maxf) maxf = f.sqnorm;

			newsols[idx] = Complex.sub(newsols[idx], f);
		}

		sols = newsols;

		if (maxf < 1e-30) break;
	}

	return sols;
}

export function solveQuarticBairstow(p, q, r, s) {
	// Solves the equation
	//  4     3     2
	// x  + px  + qx  + rx + s = 0

	let converged = false;
	let iters = 0;

	let u, v;

	while (!converged && iters < 100) {
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

export function solveQuarticAnalytic(B, C, D, E) {
	// Equation: x^4 + Bx^3 + Cx^2 + Dx + E = 0

	let a = -0.375 * B * B + C;
	let b = 0.125 * B * B * B - 0.5 * B * C + D;
	let c = -3 / 256 * B * B * B * B + 0.0625 * B * B * C - 0.25 * B * D + E;

	// Equation: u^4 + au^2 + bu + c = 0, where u = x + B / 4

	if (b === 0) {
		// biquadratic
		let uroots = solveQuadratic(a, c);
		let root0 = Complex.sqrt(uroots[0]);
		let root2 = Complex.sqrt(uroots[1]);

		return [
			Complex.add(Complex.fromReal(-B / 4), root0),
			Complex.sub(Complex.fromReal(-B / 4), root0),
			Complex.add(Complex.fromReal(-B / 4), root2),
			Complex.sub(Complex.fromReal(-B / 4), root2)
		];
	}

	// Not a biquadratic
	let y = solveCubicAnalytic(-a / 2, -c, 0.5 * (a * c - 0.25 * b * b))[0];

	let root = Complex.sqrt(Complex.sub(Complex.multReal(y, 2), Complex.fromReal(a)));
	let sum = Complex.sub(Complex.multReal(y, -2), Complex.fromReal(a));

	let u0 = Complex.multReal(
		Complex.add(
			root, 
			Complex.sqrt(
				Complex.sub(
					sum, 
					Complex.divide(
						Complex.fromReal(2 * b), 
						root
					)
				)
			)
		), 
		0.5
	);

	let u1 = Complex.multReal(
		Complex.sub(
			root, 
			Complex.sqrt(
				Complex.sub(
					sum, 
					Complex.divide(
						Complex.fromReal(2 * b), 
						root
					)
				)
			)
		), 
		0.5
	);

	let u2 = Complex.multReal(
		Complex.add(
			root, 
			Complex.sqrt(
				Complex.add(
					sum, 
					Complex.divide(
						Complex.fromReal(2 * b), 
						root
					)
				)
			)
		), 
		-0.5
	);

	let u3 = Complex.multReal(
		Complex.sub(
			root, 
			Complex.sqrt(
				Complex.add(
					sum, 
					Complex.divide(
						Complex.fromReal(2 * b), 
						root
					)
				)
			)
		), 
		-0.5
	);

	u0.a -= B / 4;
	u1.a -= B / 4;
	u2.a -= B / 4;
	u3.a -= B / 4;

	return [
		u0, 
		u1, 
		u2, 
		u3
	];
}

export function realRootsQuartic(p, q, r, s) {
	let roots = solveQuarticBairstow(p, q, r, s);
	let out = [];

	for (let i = 0; i < 4; i++) {
		if (Math.abs(roots[i].b) < 1e-10) {
			out.push(roots[i].a);
		}
	}

	out.sort();

	return out;
}