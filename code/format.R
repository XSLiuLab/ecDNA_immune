dt1 <- readxl::read_xlsx("results/Figures/tables/table1.xlsx")
myft <- flextable(dt1)
myft <- autofit(myft)
myft <- theme_box(myft)
myft
save_as_image(x = myft, path = "../tmp/test.pdf")

print(myft,preview = "pdf")

height(myft, height = .5) %>% 
  hrule(rule = "exact", part = "body")
