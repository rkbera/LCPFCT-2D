        IMPLICIT NONE
        INTEGER NPTX, NPTY
        INTEGER NXINT, NYINT
        INTEGER NPT, NWORK

        PARAMETER ( NPTX  = 500, NPTY  = 500)
        PARAMETER ( NXINT = NPTX+1, NYINT = NPTY+1)
        PARAMETER ( NPT   = 3000000, NWORK = 3000000) !change at your own risk -- RAM dependent for single processor running
