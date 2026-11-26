# Data sets

| File | Variable | n | Description | Source |
|---|---|---|---|---|
| `ncyoungmale.csv` | `proportion` | 90 | Data I. Proportion of the county population that is male and aged 15–24 years (`pctymle`) in 90 counties of North Carolina (USA) in 1987; `county` is the county identifier. | `crime4` data in the R package **wooldridge** (also `Crime` in **Ecdat**); Cornwell C, Trumbull WN. Estimating the economic model of crime with panel data. *Rev Econ Stat*. 1994;76(2):360–6. |
| `catholic.csv` | `proportion` | 46 | Data II. Proportion of Catholics (as opposed to Protestants) in 47 French-speaking provinces of Switzerland, about 1888 (`Catholic / 100`). The province of Hérens (exactly 100%) lies on the boundary of the support and is excluded. | `swiss` data set in base R (package **datasets**); Mosteller F, Tukey JW. *Data Analysis and Regression*. Reading (MA): Addison-Wesley; 1977. |

Both files can be re-created in R:

```r
# Data I
library(wooldridge); data(crime4)
ncyoungmale <- subset(crime4, year == 87, select = c(county, pctymle))
names(ncyoungmale)[2] <- "proportion"
# Data II
catholic <- swiss$Catholic / 100; catholic <- catholic[catholic < 1]
```
