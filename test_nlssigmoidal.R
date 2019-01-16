

dat <- read.table(text="Species	WP	PLC
Ep3R1	-0.44	0
Ep3R1	-1.03	0.991516
Ep3R1	-1.57	5.519779
Ep3R1	-2.11	5.535112
Ep3R1	-2.65	7.191046
Ep3R1	-3.15	16.05336
Ep3R1	-3.65	66.39068
Ep3R1	-4.14	91.567
Ep3R1	-4.65	98.54339", header=TRUE)

f <- fitplc(dat, varnames = c(PLC = "PLC", WP = "WP"), 
            model = "nls_sigmoidal")
