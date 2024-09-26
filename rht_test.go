package radiative

import (
	"fmt"
	"math"
	"testing"
)

func Pow4(T float64) float64 {
	return math.Pow(T, 4)
}

// 1D
func TestRHT1(t *testing.T) {
	const N = 20
	var (
		ak, tg [N]float64
		// input data
		tw1 = 1000.0 // K
		tw2 = 500.0  // K
		akd = 2.0
		em1 = 0.5
		em2 = 0.5
	)

	dak := akd / float64(N)
	for i := 0; i < N; i++ {
		ak[i] = (float64(i) + 0.5) * dak
	}

	tav := (tw1 + tw2) * 0.5
	for i := 0; i < N; i++ {
		tg[i] = tav
	}

	sbc := 5.6687e-8
	pai := 3.14159
	var aipo float64
	if em1 == 1.0 {
		aipo = sbc * Pow4(tw1) / pai
	}
	var ainkd float64
	if em2 == 1.0 {
		ainkd = sbc * Pow4(tw2) / pai
	}
	// 5000
L5000:
	if em1 < 1 || em2 < 1 {
		var a, b float64
		for i := 0; i < N; i++ {
			a += Pow4(tg[i]) * math.Exp(-2.0*ak[i]) * dak
			b += Pow4(tg[i]) * math.Exp(-2.0*(akd-ak[i])) * dak
		}
		a = a * 2.0 * sbc / pai
		b = b * 2.0 * sbc / pai
		if em1 < 1 {
			aipo = (em1*sbc*Pow4(tw1)/pai + (1.0-em1)*a + (1.0-em1)*
				math.Exp(-2.0*akd)*(em2*sbc*Pow4(tw2)/pai+(1.0-em2)*b)) /
				(1.0 - (1.0-em1)*(1.0-em2)*math.Pow(math.Exp(-2.0*akd), 2.0))
		}
		if em2 < 1 {
			ainkd = ((1.0-em2)*math.Exp(-2.0*akd)*(em1*sbc*Pow4(tw1)/pai+
				(1.0-em1)*a) + em2*sbc*Pow4(tw2)/pai + (1.0-em2)*b) /
				(1.0 - (1.0-em1)*(1.0-em2)*math.Pow(math.Exp(-2.0*akd), 2.0))
		}
	}
	// 49
	eps := -1.0
	for i := 0; i < N; i++ {
		tp := tg[i]
		ti := 0.0
		for j := 0; j < N; j++ {
			ti = ti + Pow4(tg[j])*math.Exp(-2.0*math.Abs(ak[j]-ak[i]))*dak
		}
		tg[i] = math.Pow(ti+pai*(aipo*math.Exp(-2.0*ak[i])+ainkd*math.Exp(-2.0*
			(akd-ak[i]))/(2.0*sbc)), 0.25)
		epsi := math.Abs(tg[i]-tp) / tg[i]
		if eps < epsi {
			eps = epsi
		}
	}
	if 1.0e-5 < eps {
		goto L5000
	}
	qw := 0.0
	for i := 0; i < N; i++ {
		qw += Pow4(tg[i]) * math.Exp(-2.0*ak[i]) * dak
	}
	qw = pai*(aipo-ainkd*math.Exp(-2.0*akd)) - 2.0*sbc*qw
	qnd := qw / (sbc * (Pow4(tw1) - Pow4(tw2)))

	for i := 0; i < N; i++ {
		fmt.Printf("%02d %.5e\n", i, tg[i])
	}
	fmt.Println(qw)
	fmt.Println(qnd)
}
