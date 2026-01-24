const phi = (1 + Math.sqrt(5)) / 2;

export function brentMinimize(f, min, max, tol, maxiters) {
	// Implementation taken from Astrogoth's Kerbal Transfer Illustrator
	const GRsq = 2 - phi; // 1 / (phi * phi)

	let a = min, b = max;
	let h = (b - a);
	let x = a * GRsq * h;
	let w = x;
	let v = w;

	let fv = f(v);
	let fw = f(w);
	let fx = f(x);

	let u = 0, fu = 0;

	let m = (a + b) / 2;

	let p = 0, q = 0, delta = 0, d = 0, e = 0;

	let it = 0;
	while (it < maxiters) {
		it++;

		if (h < tol) return [m, f(m)];

		p = (w - x) * (w - x) * (fx - fv) + (v - x) * (v - x) * (fw - fx);
		q = (w - x) * (fx - fv) + (v - x) * (fw - fx);

		delta = 0.5 * p / q;

		if (
			(q === 0) ||
			(a > u || u > b) ||
			(Math.abs(p / q) > 0.5 * Math.abs(e)) ||
			(Math.abs(e) < tol)
		) {
			e = (x < m) ? (b - x) : (a - x);
			d = GRsq * e;
		} else {
			e = d;
			d = delta;
		}

		u = x + d;
		fu = f(u);

		if (fu <= fx) {
			a = (u < x) ? a : x;
			b = (u < x) ? x : b;
			v = w;
			w = x;
			x = u;
			fv = fw;
			fw = fx;
			fx = fu;
		} else {
			a = (u < x) ? u : a;
			b = (u < x) ? b : u;

			if((fu <= fw) || (w === x)) {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else if((fu <= fv) || (v === x) || (v === w)) {
                v = u;
                fv = fu;
            }
		}

		h = b - a;
		m = (a + b) / 2;
	}

	return [m, f(m)];
}