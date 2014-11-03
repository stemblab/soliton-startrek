movie = new $blab.Movie
    aniId: "solitons"
    stackId: null
    N: 256
    h: 4e-5
    dispersion: (z) -> j*z.pow(3)

setTimeout (->
    movie.initSoliton(-1, 800)
    movie.initSoliton(1, 200)
    ), 3000
