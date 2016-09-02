

## GLM with emphasis on visualization.
## Available models include glm(...) and bayesglm(...) from `arm`
##+package. For some reason that perfect sepperation occurs, the 
##+coefficients are inflated. In this case, switch to use `bayesglm`
##+is recommended.

## ref: http://stats.stackexchange.com/questions/11109/how-to-deal-with-perfect-separation-in-logistic-regression

glmfit <-
  function(...,
           clip = c(0, 8),
           user_specified_variables,
           plotting_for_sig_only = FALSE,
           main = NULL, main.xycoord = c(0.3, 1.05),
           use.stats.glm=FALSE)
  {
  if (use.stats.glm) {
    fit <- stats::glm(...)
  } else {
    fit <- arm::bayesglm(...)
  }
  r <- arm::display(fit, digits=6, detail=TRUE)
  if (plotting_for_sig_only && all(r$p.value[2:length(r$p.value)] > 0.05))
    return(fit)
  m <- signif(exp(r$coef), 2)
  l <- signif(exp(r$coef - 2 * r$se), 2)
  u <- signif(exp(r$coef + 2 * r$se), 2)
  p <- signif(r$p.value, 2)
  if (any(is.infinite(u))) {
    warning("Inf found for coefficients, set use.stats.glm = FALSE is recommended.")
  }
  if (plotting_for_sig_only && all(l <= 1)) ## works for logit only
    return(fit)
  if (missing(user_specified_variables)) {
    variable_labels <- names(m)
  } else { ## Suppose that user_specified_variables does not include a term for intercept.
    ## I use NA for the intercept term.
    variable_labels <- c(NA, user_specified_variables)
  }
  tabletext <- list()
  tabletext[[1]] <- c("Var", "OR", "2.5% CI", "97.5% CI", "P-value")
  for (i in 2:length(m)) { ## Intercept was not included
    tabletext[[i]] <- c(variable_labels[i], m[i], l[i], u[i], p[i])
  }
  tabletext <- do.call("rbind", tabletext)
  is.summary = c(1, rep(0, 20))
  m[1] <- NA  ## Set the 1st item (i.e. intercept) for NA, empty plot for this term.
  l[1] <- NA
  u[1] <- NA
  forestplot(
    tabletext,
    m, l, u,
    clip = clip,
    zero = 1,
    is.summary = is.summary,
    boxsize = 0.5,
    align = rep('l',5),
    xlab = "OR (95% CI)",
    col = meta.colors(box = "black", line = "black", summary = "black", zero = "black")
  )
  if (!is.null(main)) {
    x <- main.xycoord[1]
    y <- main.xycoord[2]
    text(x, y, main, xpd = TRUE)
  }
  invisible(fit)
}

library(arm)
library(rmeta)

## An example of glmfit in regressing SMG mutation status versus mutational exposures
stad_SMGs_vs_MSig <- function() {
  setwd(
    "/Users/lixiangchun/Public/WorkSpace/Project/DigestiveSystemCancer/Analysis/iCGA_PanelOfNormals/StomachPublishedDataOnlyReanalyses/mutated_genes_vs_mutational_signatures"
  )
  load('SMG_vs_MSig.RData')
  pdf(
    "RegularMutated_Smgs_vs_mutational_signatures.pdf",
    width = 7,
    height = 4
  )
  for (gene in regular_mutated_smgs) {
    f <- sprintf("%s ~ msig1 + msig2 + msig3 + msig4 + msig5 + msig6 + msig7", gene)
    fit = glmfit(
      as.formula(f),
      data = regular,
      family = binomial,
      clip = c(0, 9),
      main = sprintf("%s vs. mutational exposures", gene),
      use.stats.glm = FALSE, plotting_for_sig_only = TRUE
    )
  }
  dev.off()
  pdf(
    "HyperMutated_Smgs_vs_mutational_signatures.pdf",
    width = 7,
    height = 4
  )
  for (gene in hyper_mutated_smgs) {
    f <- sprintf("%s ~ msig1 + msig2 + msig3 + msig4 + msig5 + msig6 + msig7", gene)
    fit = glmfit(
      as.formula(f),
      data = hyper,
      family = binomial,
      clip = c(0, 9),
      main = sprintf("%s vs. mutational exposures", gene),
      use.stats.glm = FALSE, plotting_for_sig_only = TRUE
    )
  }
  dev.off()
}

##stad_SMGs_vs_MSig()
