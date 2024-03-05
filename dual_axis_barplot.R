library(ggplot2)

# create test data
crispra <- sample(1:100, 50)
crispra <- crispra[order(crispra)]
crisprn <- sample(1:100, 25)
crisprn <- crisprn[order(-crisprn)]
plotdat <- data.frame(x = 1:length(crispra),
                      crispra = crispra,
                      crisprn = c(crisprn, rep(NA, 25)))

# plot
ggplot(plotdat, aes(x = x)) +
  geom_bar(aes(y = crispra), stat = "identity", fill = "firebrick") + # bars for first crispr set
  geom_segment(aes(xend = x, y = 100, yend = 100-crisprn), 
               size = 5, stat = "identity", col = "forestgreen") + # bars for second - have to use geom_segment to get custom start/end points
  # note you need to change the size paramter to match between bars & segments depending on size of plotting window
  scale_y_continuous(sec.axis = sec_axis(~. - 100, name = "CRISPRn"), name = "CRISPa") # add second axis 
