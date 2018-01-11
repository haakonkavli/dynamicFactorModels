dfm <- function(dat_input = NULL,
                n_factors = 1,
                model_structure = c("diagonal and equal")){


    # Implements MARSS dfm

    ## Inputs:
    # dat_input -- X table with dimension (Nx variables, m Examples)
    # n_factors -- Number of Y variables to extract, Ny
    # models_structure -- Restrictions on covariance matrix, see MARSS help

        if (is.null(dat_input)) {
                stop('Must pass in data in dat_input..')
        }

        n_variables = dim(dat_input)[1]
        n_obs = dim(dat_input)[2]

        # Prepare output list:
        output = list(data = data.table(),
                      fit = data.table(),
                      trends = data.table(),
                      loads = data.table())

        # Standardise data:
        sigma = sqrt(apply(dat_input, 1, var, na.rm = TRUE))
        y_bar = apply(dat_input, 1, mean, na.rm = TRUE)
        z_scores = (dat_input - y_bar) * (1/sigma) # Z scored data
        rownames(z_scores) = rownames(dat_input)

        # Set control params
        control_pars = list(minit = 500,
                            maxit = 10000,
                            allow.degen = FALSE)

        # set up forms of R matrices
        dfm = data.table()

        print(paste(Sys.time(),": Estimating DFM"))

        dfm_pars = list(A = "zero",
                        R = model_structure,
                        m = n_factors)

        dfm_estimate = MARSS(z_scores,
                             model = dfm_pars,
                             control = control_pars,
                             form = "dfa",
                             z.score = TRUE)

        dfm = rbind(dfm,
                    data.table(`Model Structure` = model_structure,
                               `Factors` = n_factors,
                               `Log Likelihood` = dfm_estimate$logLik,
                               `Parameters` = dfm_estimate$num.params,
                               `AIC` = dfm_estimate$AICc,
                               stringsAsFactors = FALSE))

        # Extract rotated trends (factors)
        if (n_factors > 1) {

                H_inverse = varimax(coef(dfm_estimate,
                                         type = "matrix")$Z)$rotmat

                # rotate factor loadings
                rotated_parameteres = coef(dfm_estimate,
                                           type = "matrix")$Z %*% H_inverse

                # rotate trends
                trends = solve(H_inverse) %*% dfm_estimate$states

        } else {

                rotated_parameteres <- coef(dfm_estimate, type = "matrix")$Z
                trends <- dfm_estimate$states

        }

        # Prepare output:
        ts_trends = t(trends)
        estimated_parameters = coef(dfm_estimate, type = "matrix")
        fitted_data = estimated_parameters$Z %*% dfm_estimate$states +
                matrix(estimated_parameters$A,
                       nrow = n_variables,
                       ncol = n_obs)

        rownames(fitted_data) <- rownames(z_scores)

        colnames(rotated_parameteres) <- colnames(ts_trends)

        rownames(rotated_parameteres) <- rownames(t(fitted_data))

        output$input_z_scores <- z_scores
        output$fitted_data <- t(fitted_data)
        output$trends <- ts_trends
        output$loads <- rotated_parameteres
        output$sigma_mean <- rbind(sigma,y_bar)
        return(output)
}
