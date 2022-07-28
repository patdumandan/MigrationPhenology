require(dplyr)
require(lubridate)
require(tidyr)

gos=GosMts%>%select(-TOTAL)%>%
  pivot_longer(cols=3:27, names_to = "Species", values_to="Count")%>%
  mutate(year=year(Date), month=month(Date), date=day(Date))%>%
  filter(!is.na(Count))

gos$Julian=as.POSIXlt(gos$Date)$yday

gos=gos%>%group_by(year)%>%mutate(yr_tot=cumsum(Count))
