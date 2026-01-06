import { Vector3 } from './vector.js';

function hypergeometric(x, tol) {
	let y = 1;
	let term = 1;

	for (let j = 0; j < 100; j++) {
		term *= (j + 3) / (j + 2.5) * x;

		if (Math.abs(term) < tol) break;

		y += term;
	}

	return y;
}

function x2tofBattin(x, N, lambda) {
	let K = lambda * lambda;
	let E = x * x - 1;
	let rho = Math.abs(E);
	let z = Math.sqrt(1 + K * E);

	let eta = z - lambda * x;
	let S1 = 0.5 * (1 - lambda - x * eta);
	let Q = 4 / 3 * hypergeometric(S1, 1e-11);

	return (eta * eta * eta * Q + 4 * lambda * eta) / 2 + N * Math.PI / Math.pow(rho, 1.5);
}

function x2tofLagrange(x, N, lambda) {
	let a = 1 / (1 - x * x);

	let tof;
	if (a > 0) {
		let alpha = 2 * Math.acos(x);
		let beta = 2 * Math.asin(Math.sqrt(lambda * lambda / a));

		if (lambda < 0) beta = -beta;

		tof = a * Math.sqrt(a) * ((alpha - Math.sin(alpha)) - (beta - Math.sin(beta)) + 2 * Math.PI * N) / 2;
	} else {
		let alpha = 2 * Math.acosh(x);
		let beta = 2 * Math.asinh(Math.sqrt(-lambda * lambda / a));

		if (lambda < 0) beta = -beta;

		tof = a * Math.sqrt(-a) * ((alpha - Math.sinh(alpha)) - (beta - Math.sinh(beta))) / 2;
	}

	return tof;
}

function x2tofLancaster(x, N, lambda) {
	let K = lambda * lambda;
	let E = x * x - 1;
	let rho = Math.abs(E);
	let z = Math.sqrt(1 + K * E);

	let y = Math.sqrt(rho);
	let g = x * z - lambda * E;
	let d;

	if (E < 0) {
		let l = Math.acos(g);
		d = N * Math.PI + l;
	} else {
		let f = y * (z - lambda * x);
		d = Math.log(f + g);
	}

	return (x - lambda * z - d / y) / E;
}

function x2tof(x, N, lambda) {
	const battin = 0.01;
	const lagrange = 0.2;

	let dist = Math.abs(x - 1);

	if (dist < battin) return x2tofBattin(x, N, lambda);
	if (dist < lagrange) return x2tofLagrange(x, N, lambda);
	return x2tofLancaster(x, N, lambda);
}

function dTdx(x, T, lambda) {
	let l2 = lambda * lambda;
	let l3 = lambda * l2;
	let umx2 = 1 - x * x;
	let y = Math.sqrt(1 - l2 * umx2);
	let y2 = 1 - l2 * umx2;
	let y3 = y2 * y;

	let DT = (3 * T * x - 2 + 2 * l3 * x / y) / umx2;
	let DDT = (3 * T + 5 * x * DT + 2 * (1 - l2) * l3 / y3) / umx2;
	let DDDT = (7 * x * DDT + 8 * DT - 6 * (1 - l2) * l2 * l3 * x / y3 / y2) / umx2;

	return [DT, DDT, DDDT];
}

function householder(x0, T, N, lambda, eps, maxiters) {
	let x = x0;

	for (let i = 0; i < maxiters; i++) {
		let tof = x2tof(x, N, lambda);

		let derivatives = dTdx(x, tof, lambda);
		let DT = derivatives[0];
		let DDT = derivatives[1];
		let DDDT = derivatives[2];

		let delta = tof - T;
		let DT2 = DT * DT;

		let xnew = x - delta * (DT2 - delta * DDT / 2) / (DT * (DT2 - delta * DDT) + DDDT * delta * delta / 6);
		let err = Math.abs(x - xnew);

		if (err < eps) break;

		x = xnew;
	}

	return x;
}

export function lambert(r1, r2, tof, mu, maxrevs, retrograde) {
	// An implementation of Izzo's algorithm based off of pykep's implementation
	let chord = Vector3.sub(r2, r1).norm;
	let dist1 = r1.norm;
	let dist2 = r2.norm;

	let s = 0.5 * (chord + dist1 + dist2);

	let ir1 = Vector3.normalize(r1);
	let ir2 = Vector3.normalize(r2);

	let ih = Vector3.cross(ir1, ir2);

	let lambda2 = 1 - chord / s;
	let lambda = Math.sqrt(lambda2);

	let it1, it2;

	if ((ih.z < 0) !== retrograde) {
		lambda *= -1;

		it1 = Vector3.normalize(Vector3.cross(ir1, ih));
		it2 = Vector3.normalize(Vector3.cross(ir2, ih));
	} else {
		it1 = Vector3.normalize(Vector3.cross(ih, ir1));
		it2 = Vector3.normalize(Vector3.cross(ih, ir2));
	}

	let T = Math.sqrt(2 * mu / s) * tof / s;
	let lambda3 = lambda * lambda2;

	//console.log(chord, s, lambda, T);

	// Find x based on T, lambda
	let Nmax = Math.trunc(T / Math.PI);
	let T00 = Math.acos(lambda) + lambda * Math.sqrt(1 - lambda2);
	let T0 = (T00 + maxrevs * Math.PI);
	let T1 = 2 / 3 * (1 - lambda3);
	
	if ((Nmax > 0) && (T < T0)) {
		let err = 1;
		let Tmin = T0;
		let xold = 0;
		let xnew = 0;

		for (let i = 0; i <= 12; i++) {
			let derivatives = dTdx(xold, Tmin, lambda);

			let DT = derivatives[0];
			let DDT = derivatives[1];
			let DDDT = derivatives[2];

			if (DT !== 0) {
				xnew = xold - DT * DDT / (DDT * DDT - DT * DDDT / 2);
			}

			err = Math.abs(xold - xnew);

			if (err < 1e-13) break;

			Tmin = x2tof(xnew, Nmax, lambda);
			xold = xnew;
		}
		if (Tmin > T) {
			Nmax--;
		}
	}
	if (Nmax > maxrevs) {
		Nmax = maxrevs;
	}

	let xs = Array(2 * Nmax + 1);

	if (T >= T00) {
		xs[0] = (T00 - T) / (T - T00 + 4);
	} else if (T <= T1) {
		xs[0] = T1 * (T1 - T) / (2 / 5 * (1 - lambda2 * lambda3) * T) + 1;
	} else {
		xs[0] = Math.pow(T / T00, Math.log(2) / Math.log(T1 / T00)) - 1;
	}

	xs[0] = householder(xs[0], T, 0, lambda, 1e-5, 15);

	for (let i = 1; i <= Nmax; i++) {
		let tmp = Math.pow((i * Math.PI + Math.PI) / (8 * T), 2 / 3);
		xs[2 * i - 1] = (tmp - 1) / (tmp + 1);
		xs[2 * i - 1] = householder(xs[2 * i - 1], T, i, lambda, 1e-8, 15);

		tmp = Math.pow((8 * T) / (i * Math.PI), 2 / 3);
		xs[2 * i] = (tmp - 1) / (tmp + 1);
		xs[2 * i] = householder(xs[2 * i], T, i, lambda, 1e-8, 15);
	}

	//console.log(xs);

	// Get v1, v2
	let sols = [];

	let gamma = Math.sqrt(mu * s / 2);
	let rho = (dist1 - dist2) / chord;
	let sigma = Math.sqrt(1 - rho * rho);
	for (let i = 0; i < xs.length; i++) {
		let y = Math.sqrt(1 - lambda2 + lambda2 * xs[i] * xs[i]);

		let vr1 = gamma * ((lambda * y - xs[i]) - rho * (lambda * y + xs[i])) / dist1;
		let vr2 = -gamma * ((lambda * y - xs[i]) + rho * (lambda * y + xs[i])) / dist2;

		let vt = gamma * sigma * (y + lambda * xs[i]);

		let vt1 = vt / dist1;
		let vt2 = vt / dist2;

		sols.push([
			Vector3.add(
				Vector3.mult(ir1, vr1),
				Vector3.mult(it1, vt1)
			), 
			Vector3.add(
				Vector3.mult(ir2, vr2),
				Vector3.mult(it2, vt2)
			)
		]);
	}

	return sols;
}