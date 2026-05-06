par(mfrow=c(1,2))

hist(crop$PercSand, breaks = 20)
hist(crop$PercSand %>% log(base = 10), breaks = 20)

hist(crop$PercTotNitro, breaks = 20)
hist(crop$PercTotNitro %>% log(base = 10), breaks = 20)

hist(crop$PercSOC, breaks = 20)
hist(crop$PercSOC %>% log(base = 10), breaks = 20)

hist(crop$PercTotPhos, breaks = 20)
hist(crop$PercTotPhos %>% log(base = 10), breaks = 20)

hist(crop$PercClay, breaks = 20)
hist(crop$PercClay %>% log(base = 10), breaks = 20)

hist(crop$PercSilt, breaks = 20)
hist(crop$PercSilt %>% log(base = 10), breaks = 20)

hist(crop$PH, breaks = 20)
hist(crop$PH %>% log(base = 10), breaks = 20)


# Bimodalità di BulkDensity
plot_x_vs_y(crop, 'BulkDensity')





