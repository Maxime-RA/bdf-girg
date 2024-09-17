from girg_sampling import girgs


def test_generate_girg():
    alpha = 100
    dim = 4
    n = 135
    deg = 4.2

    w = girgs.generateWeights(n, 1.5, seed=41)
    assert len(w) == n

    p = girgs.generatePositions(n, dim, seed=42)
    assert len(p) == n
    assert len(p[0]) == dim

    ws = girgs.scaleWeights(w, deg, dim, alpha)
    w = [x * ws for x in w]

    e = girgs.generateEdges(w, p, alpha, seed=43)
    assert len(e) > 0.8 * (n * deg / 2)
    assert len(e) < 1.2 * (n * deg / 2)


