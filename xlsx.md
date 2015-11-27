
# Writing files with xlsx package

## The proper way

I will put the shortcut later in the document but the
full way to write excel files from R follows.

First you have to create a workbook object

`wb <- createWorkbook(type="xlsx")`

And then create sheets to add to it

`sheet <- createSheet(wb, sheetName = "testsheet")`

You can then add a data frame to the sheet

`d <- data.frame(f=gl(2,100,labels = c("A","B")), x=rnorm(100)))`

Add the newly created data frame

`addDataFrame(d, sheet)`

You can add another sheet

```
d$y <- rbinom(100,10,0.5)
sheet <- createSheet(wb, sheetName = "sheet 2")
addDataFrame(d, sheet)
```

Finally write the sheet to a file


`saveWorkbook(wb, file="Overlaps.xlsx")`

## The Quick Way

`write.xlsx(d, "out.xlsx", sheetName="Sheet 1", col.names=TRUE, row.names=FALSE, append=FALSE)`

