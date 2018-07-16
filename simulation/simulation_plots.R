
site_range = c(10, 100, 1000, 10000, 100000)

n = c(0.1723604659288071, 0.5092162108359479, 0.9192619172938433, 0.9913056380248052, 0.9991131920452341)
qp = c(0.13279807209978292, 0.49568256543353106, 0.9205698535682051, 0.991372434456825, 1)

pdf("~/Desktop/bla.pdf") 
plot(site_range, n, log="x", pch=19, ylim=c(0,1), col="blue")
lines(site_range, n, col="blue")
points(site_range, qp, pch=19, col="red")
lines(site_range, qp, col="red")
dev.off()

n = c(0.9965352855243402, 0.7478084999261679, 0.10356482884972121)
qp = c(0.9967598240044038, 0.7422261301738312, 0.1)
tissue_range = c(10, 100, 1000)
plot(tissue_range, n, log="x", pch=19, ylim=c(0,1), col="blue")
lines(tissue_range, n, col="blue")
points(tissue_range, qp, pch=19, col="red")
lines(tissue_range, qp, col="red")

n = c(0.9181521468164512, 0.9905975162556996, 0.9991264930706049)
qp = c(0.9240631181496793, 0.9908250548647307, 1)
depth_range = c(10, 100, 1000)
plot(depth_range, n, log="x", pch=19, ylim=c(0,1), col="blue")
lines(depth_range, n, col="blue")
points(depth_range, qp, pch=19, col="red")
lines(depth_range, qp, col="red")
grid(nx = 7, ny = nx, col = "lightgray", lty = "dotted",
lwd = par("lwd"), equilogs = TRUE)



[0.9908998773369697, 0.981560137535239, 0.9582227658700151, 0.5885503623735578]
[0.9906951293416658, 0.9822650570487778, 0.9573988193139126, 0.5769283228569685]





