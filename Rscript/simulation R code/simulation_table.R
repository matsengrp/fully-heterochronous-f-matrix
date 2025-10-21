#source("diag_top_down_final.R")
#source("bernoulli_splitting_final.R")
#source("sampling_coales_final.R")


library(xtable)

## --- numbers (fill/modify as needed) ---
coal_T5  <- c(coales_sample_1000_5$avg_cherry, coales_sample_1000_5$avg_int_length, coales_sample_1000_5$avg_total_length)
coal_T20 <- c(coales_sample_1000_20$avg_cherry, coales_sample_1000_20$avg_int_length, coales_sample_1000_20$avg_total_length)
coal_T50 <- c(coales_sample_1000_50$avg_cherry, coales_sample_1000_50$avg_int_length, coales_sample_1000_50$avg_total_length)

diag_T5  <- c(diag_sample_1000_5$avg_cherry, diag_sample_1000_5$avg_int_length, diag_sample_1000_5$avg_total_length)
diag_T20 <- c(diag_sample_1000_20$avg_cherry, diag_sample_1000_20$avg_int_length, diag_sample_1000_20$avg_total_length)
diag_T50 <- c(diag_sample_1000_50$avg_cherry, diag_sample_1000_50$avg_int_length, diag_sample_1000_50$avg_total_length)


statlab <- c( "$N_C$", "$\\bar{ L^{I}}$", "$\\bar{ L^{T}}$")

df <- data.frame(
  `  ` = statlab,
  coal_T5, coal_T20, coal_T50,
  diag_T5, diag_T20, diag_T50,
  check.names = FALSE
)

## --- no separators: simple alignment (ncol(df)=8 => length 9) ---
xt <- xtable(
  df, caption = " ", label = "tab:sim.coales_diag",
  align = c("l","l","r","r","r","r","r","r"),
  digits = c( 0,0, rep(2, 6))
)

## --- two-row header WITHOUT rules/lines ---
hdr <- paste0(
  "  & \\multicolumn{3}{c}{coalescence}",
  " & \\multicolumn{3}{c}{diagonal top-down} \\\\[-2pt]\n",
  "  $n=100$ & $|T|=5$ & $|T|=20$ & $|T|=50$",
  " & $|T|=5$ & $|T|=20$ & $|T|=50$ \\\\ \n"
)

add <- list(
  pos = list(-1),             # insert only before the first row
  command = hdr
)

print(xt,
      include.rownames = FALSE,
      include.colnames = FALSE,
      sanitize.text.function = identity,
      add.to.row = add,
      floating = FALSE)


## --- numbers (fill/modify as needed) ---
B_a10_b1_T5  <- c(bernoulli_1000_5.10_1$avg_cherry, bernoulli_1000_5.10_1$avg_int_length, bernoulli_1000_5.10_1$avg_total_length)
B_a10_b1_T20 <- c(bernoulli_1000_20.10_1$avg_cherry, bernoulli_1000_20.10_1$avg_int_length, bernoulli_1000_20.10_1$avg_total_length)
B_a10_b1_T50 <- c(bernoulli_1000_50.10_1$avg_cherry, bernoulli_1000_50.10_1$avg_int_length, bernoulli_1000_50.10_1$avg_total_length)

B_a10_b10_T5  <- c(bernoulli_1000_5.10_10$avg_cherry, bernoulli_1000_5.10_10$avg_int_length, bernoulli_1000_5.10_10$avg_total_length)
B_a10_b10_T20 <- c(bernoulli_1000_20.10_10$avg_cherry, bernoulli_1000_20.10_10$avg_int_length, bernoulli_1000_20.10_10$avg_total_length)
B_a10_b10_T50 <- c(bernoulli_1000_50.10_10$avg_cherry, bernoulli_1000_50.10_10$avg_int_length, bernoulli_1000_50.10_10$avg_total_length)

B_a1_b10_T5  <- c(bernoulli_1000_5.1_10$avg_cherry, bernoulli_1000_5.1_10$avg_int_length, bernoulli_1000_5.1_10$avg_total_length)
B_a1_b10_T20 <- c(bernoulli_1000_20.1_10$avg_cherry, bernoulli_1000_20.1_10$avg_int_length, bernoulli_1000_20.1_10$avg_total_length)
B_a1_b10_T50 <- c(bernoulli_1000_50.1_10$avg_cherry, bernoulli_1000_50.1_10$avg_int_length, bernoulli_1000_50.1_10$avg_total_length)


df_2<- data.frame(
  `  ` = statlab,
  B_a10_b1_T5, B_a10_b1_T20, B_a10_b1_T50,
  B_a10_b10_T5, B_a10_b10_T20, B_a10_b10_T50,
  B_a1_b10_T5,B_a1_b10_T20,B_a1_b10_T50,
  check.names = FALSE
)


xt2 <- xtable(
  df_2, caption = " ", label = "tab:sim.coales_diag",
  align = c("l","l","r","r","r","r","r","r","r","r","r"),
  digits = c( 0,0, rep(2, 9))
)

hdr2 <- paste0(
  "  & \\multicolumn{3}{c}{$\\alpha=10,\\beta=1$}",
  " & \\multicolumn{3}{c}{$\\alpha=10,\\beta=1$}& ",
  "\\multicolumn{3}{c}{$\\alpha=10,\\beta=1$} \\\\ \n",
  "  $n=100$ & $|T|=5$ & $|T|=20$ & $|T|=50$",
  " & $|T|=5$ & $|T|=20$ & $|T|=50$", 
  " & $|T|=5$ & $|T|=20$ & $|T|=50$ \\\\ \n"
)

add2 <- list(
  pos = list(-1),             # insert only before the first row
  command = hdr2
)

print(xt2,
      include.rownames = FALSE,
      include.colnames = FALSE,
      sanitize.text.function = identity,
      add.to.row = add2,
      floating = FALSE)
