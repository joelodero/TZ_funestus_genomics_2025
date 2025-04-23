library(ggplot2)
library(data.table)

west <- fread('/Users/dennistpw/Projects/funestus_tz/stairwayplot/west_rift/west_rift.final.summary')
east <- fread('/Users/dennistpw/Projects/funestus_tz/stairwayplot/east_rift.final.summary')

west$pop <- 'Rift West'
east$pop <- 'Rift East'

both <- rbind(west, east)

ggplot(both, aes(x=log10(year),y=log10(Ne_median), colour=pop))+
  geom_line(size=1)+
  geom_ribbon(aes(ymin = log10(`Ne_2.5%`), ymax = log10(`Ne_97.5%`), fill=pop), alpha = 0.1) +
  theme_classic()+
  xlim(log10(100), log10(1000000000))+
  labs(x='Years before present (log10)', y='Median Ne (log10)', fill = 'Population')
