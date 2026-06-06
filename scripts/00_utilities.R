# =============================================================================
# 00_utilities.R  —  Funzioni condivise tra gli script del progetto
# =============================================================================
#
# Struttura:
#   1. Helper di progetto  (setup_rtools, carica_dati, fit_field)
#   2. EDA generica        (analizza_dataset, correlazione_*, grafico_*, plot_*)
#   3. Visualizzazioni     (plot_depth_profiles, plot_texture_triangle)
# =============================================================================


# ── 1. HELPER DI PROGETTO ─────────────────────────────────────────────────────

# Configura PATH per Rtools su Windows (necessario per cmdstanr).
setup_rtools <- function() {
  if (.Platform$OS.type != "windows") return(invisible(NULL))
  rtools_path <- Sys.getenv("RTOOLS45_HOME", unset = "C:/rtools45")
  Sys.setenv(
    PATH = paste(file.path(rtools_path, "ucrt64/bin"),
                 file.path(rtools_path, "usr/bin"),
                 Sys.getenv("PATH"), sep = ";"),
    RTOOLS44_HOME = rtools_path,
    RTOOLS45_HOME = rtools_path
  )
  invisible(NULL)
}

# Salva una figura in output/figures/ come PDF.
# Usata da tutti gli script figura (03, 12, 18, 19) in luogo della definizione locale.
save_fig <- function(fname, p, w = 16, h = 9, u = "cm") {
  ggplot2::ggsave(file.path(here::here("output", "figures"), fname),
                 plot = p, width = w, height = h, units = u, device = "pdf")
  cat(sprintf("  [fig] Salvato: %s\n", fname))
}

# Carica dati.rds e applica le trasformazioni standard (log-risposte, scaling
# di logBottom/Texture1/Texture2/BulkDensity/PH, Field come factor).
# Usata da tutti gli script modello (07-22) in luogo del blocco dati ripetuto.
carica_dati <- function() {
  readRDS(here::here("data", "dati.rds")) |>
    dplyr::mutate(dplyr::across(c(OnFarm, Irrigate, Fertilised, N_Natural),
                                ~ as.integer(as.character(.x)))) |>
    dplyr::mutate(
      logSOC    = log(PercSOC),
      logN      = log(PercTotNitro),
      logP      = log(PercTotPhos),
      logBottom = log(Bottom)
    ) |>
    dplyr::mutate(dplyr::across(c(logBottom, Texture1, Texture2, BulkDensity, PH),
                                ~ c(scale(.x)))) |>
    dplyr::mutate(Field = factor(Field))
}

# OLS per campo: log_y_var ~ logBottom centrato.
# Restituisce tibble con Field, n_obs, int (intercetta), slope, r2, slope_se.
fit_field <- function(df, y_var) {
  mu_logB <- mean(log(df$Bottom))
  if (!"logBottom_c" %in% names(df))
    df <- dplyr::mutate(df, logBottom_c = log(Bottom) - mu_logB)
  df |>
    dplyr::group_by(Field) |>
    dplyr::summarise(
      n_obs    = dplyr::n(),
      fit      = list(lm(stats::reformulate("logBottom_c", response = y_var),
                         data = dplyr::cur_data())),
      .groups  = "drop"
    ) |>
    dplyr::mutate(
      int      = purrr::map_dbl(fit, ~ stats::coef(.x)[1]),
      slope    = purrr::map_dbl(fit, ~ stats::coef(.x)[2]),
      r2       = purrr::map_dbl(fit, ~ summary(.x)$r.squared),
      slope_se = purrr::map_dbl(fit, ~ summary(.x)$coefficients[2, 2])
    ) |>
    dplyr::select(-fit)
}


# ── 2. EDA GENERICA ───────────────────────────────────────────────────────────

# Funzioni

analizza_dataset <- function(data) {
  
  # Validazione input più flessibile
  if (!is.data.frame(data) && !inherits(data, c("tbl_df", "data.table"))) {
    stop("L'input deve essere un data.frame, tibble o data.table")
  }
  
  # Gestione dataset vuoto
  if (ncol(data) == 0) {
    return(data.frame(
      variabile = character(0),
      tipo = character(0),
      campo_variazione = character(0),
      num_na = integer(0),
      num_unique = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # Pre-alloca il data.frame risultato per efficienza
  n_cols <- ncol(data)
  risultato <- data.frame(
    variabile = character(n_cols),
    tipo = character(n_cols),
    campo_variazione = character(n_cols),
    num_na = integer(n_cols),
    num_unique = integer(n_cols),
    stringsAsFactors = FALSE
  )
  
  # Funzione helper per campo di variazione
  calcola_campo_variazione <- function(col_data, num_unique, num_na) {
    # Se tutti i valori sono NA
    if (num_na == length(col_data)) {
      return("Tutti NA")
    }
    
    if (is.numeric(col_data)) {
      min_val <- min(col_data, na.rm = TRUE)
      max_val <- max(col_data, na.rm = TRUE)
      return(paste0("[", min_val, " - ", max_val, "]"))
      
    } else if (is.factor(col_data) || is.character(col_data)) {
      return(paste0(num_unique, " categorie"))
      
    } else if (is.logical(col_data)) {
      valori_unici <- unique(col_data[!is.na(col_data)])
      return(paste(valori_unici, collapse = ", "))
      
    } else if (inherits(col_data, "Date")) {
      min_date <- min(col_data, na.rm = TRUE)
      max_date <- max(col_data, na.rm = TRUE)
      # Gestione sicura per date
      if (is.na(min_date) || is.na(max_date)) {
        return("Tutti NA")
      } else {
        return(paste0("[", min_date, " - ", max_date, "]"))
      }
      
    } else {
      return(paste0(num_unique, " valori unici"))
    }
  }
  
  # Itera su ogni colonna (usando seq_len per sicurezza)
  for (i in seq_len(ncol(data))) {
    col_name <- names(data)[i]
    col_data <- data[[i]]
    
    # Calcoli una volta sola
    tipo <- class(col_data)[1]
    num_na <- sum(is.na(col_data))
    num_unique <- length(unique(col_data[!is.na(col_data)]))
    campo_var <- calcola_campo_variazione(col_data, num_unique, num_na)
    
    # Assegnazione diretta invece di rbind
    risultato[i, ] <- list(
      variabile = col_name,
      tipo = tipo,
      campo_variazione = campo_var,
      num_na = num_na,
      num_unique = num_unique
    )
  }
  
  return(risultato)
}

correlazione_x_y <- function(data, y_col = 'y', metodo = "pearson", 
                             plot = T, titolo = "Correlazioni X vs Y") {
  
  # Validazione input - accetta data.frame, tibble, data.table
  if (!is.data.frame(data)) {
    stop("L'input deve essere un data.frame, tibble o data.table")
  }
  
  if (!y_col %in% names(data)) {
    stop("La variabile Y specificata non esiste nel dataset")
  }
  
  if (!is.numeric(data[[y_col]])) {
    stop("La variabile Y deve essere numerica")
  }
  
  # Validazione metodo
  metodi_validi <- c("pearson", "spearman", "kendall")
  if (!metodo %in% metodi_validi) {
    stop(paste("Metodo deve essere uno tra:", paste(metodi_validi, collapse = ", ")))
  }
  
  # Verifica disponibilità ggplot2 se necessario
  if (plot && !requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 non disponibile. Il plot non verrà generato.")
    plot <- FALSE
  }
  
  # Seleziona solo variabili numeriche (escludendo Y)
  var_numeriche <- sapply(data, is.numeric)
  x_vars <- names(data)[var_numeriche & names(data) != y_col]
  
  if (length(x_vars) == 0) {
    stop("Non ci sono variabili X numeriche nel dataset")
  }
  
  # Pre-alloca dataframe per efficienza
  n_vars <- length(x_vars)
  correlazioni <- data.frame(
    variabile_x = x_vars,
    correlazione = numeric(n_vars),
    p_value = numeric(n_vars),
    significativa = logical(n_vars),
    stringsAsFactors = FALSE
  )
  
  y_data <- data[[y_col]]
  
  # Calcola correlazioni - riempie vettori pre-allocati
  for (i in seq_along(x_vars)) {
    x_var <- x_vars[i]
    x_data <- data[[x_var]]
    
    # Rimuovi coppie con NA
    complete_cases <- complete.cases(x_data, y_data)
    n_complete <- sum(complete_cases)
    
    if (n_complete < 3) {
      # Non abbastanza osservazioni
      correlazioni$correlazione[i] <- NA
      correlazioni$p_value[i] <- NA
      correlazioni$significativa[i] <- NA
    } else {
      # Test di correlazione
      tryCatch({
        test_result <- cor.test(x_data[complete_cases], y_data[complete_cases], 
                                method = metodo)
        
        correlazioni$correlazione[i] <- test_result$estimate
        correlazioni$p_value[i] <- test_result$p.value
        correlazioni$significativa[i] <- test_result$p.value < 0.05
      }, error = function(e) {
        correlazioni$correlazione[i] <<- NA
        correlazioni$p_value[i] <<- NA
        correlazioni$significativa[i] <<- NA
        warning(paste("Errore nel calcolo della correlazione per", x_var, ":", e$message))
      })
    }
  }
  
  # Ordina per correlazione assoluta decrescente (gestendo NA)
  correlazioni <- correlazioni[order(abs(correlazioni$correlazione), decreasing = TRUE, na.last = TRUE), ]
  rownames(correlazioni) <- NULL
  
  # Crea e stampa plot se richiesto
  if (plot && nrow(correlazioni) > 0) {
    # Rimuovi NA per il plot
    cor_plot_data <- correlazioni[!is.na(correlazioni$correlazione), ]
    
    if (nrow(cor_plot_data) > 0) {
      p <- ggplot2::ggplot(cor_plot_data, 
                           ggplot2::aes(x = reorder(variabile_x, abs(correlazione)), 
                                        y = correlazione, 
                                        fill = significativa)) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "lightgray"),
                                   name = "Significativa\n(p < 0.05)",
                                   na.value = "darkgray") +
        ggplot2::ylim(-1, 1) +
        ggplot2::coord_flip() +
        ggplot2::labs(title = titolo,
                      subtitle = paste("Metodo:", metodo),
                      x = "Variabili X",
                      y = "Correlazione") +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5))
      
      print(p)
    }
  }
  
  return(correlazioni)
}

plot_corrplot <- function(cor_matrix, titolo, metodo) {
  if (requireNamespace("corrplot", quietly = TRUE)) {
    tryCatch({
      corrplot::corrplot(cor_matrix, 
                         method = "color",
                         type = "upper",
                         order = "original",  # Evita problemi con hclust
                         tl.col = "black",
                         tl.srt = 45,
                         diag = FALSE,
                         title = paste(titolo, "-", metodo),
                         mar = c(0,0,2,0),
                         # Rimuovi coefficienti per plot più puliti
                         addCoef.col = if(ncol(cor_matrix) <= 10) "black" else NULL)
    }, error = function(e) {
      warning("Errore nel corrplot, uso imageplot base")
      plot_imageplot_matrix(cor_matrix, titolo, metodo)
    })
  } else {
    warning("Package corrplot non disponibile, uso imageplot base")
    plot_imageplot_matrix(cor_matrix, titolo, metodo)
  }
}
plot_imageplot_matrix <- function(cor_matrix, titolo, metodo) {
  # Salva parametri grafici originali
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  p <- ncol(cor_matrix)
  var_names <- colnames(cor_matrix)
  
  # Calcola dimensioni appropriate
  cex_size <- max(0.3, min(1, 15/p))
  
  # Imposta margini in base al numero di variabili
  bottom_margin <- max(3, min(15, max(nchar(var_names)) * 0.3))
  left_margin <- max(3, min(15, max(nchar(var_names)) * 0.3))
  
  par(mar = c(bottom_margin, left_margin, 4, 4))
  
  # Crea imageplot
  image(1:p, 1:p, cor_matrix,
        col = colorRampPalette(c("red", "white", "blue"))(100),
        zlim = c(-1, 1),
        xlab = "", ylab = "",
        main = paste(titolo, "\n(Metodo:", metodo, ")"),
        axes = FALSE)
  
  # Aggiungi etichette solo se leggibili
  if (p <= 50) {
    axis(1, at = 1:p, labels = var_names, las = 2, cex.axis = cex_size)
    axis(2, at = 1:p, labels = var_names, las = 1, cex.axis = cex_size)
  } else {
    # Per matrici molto grandi, mostra solo alcuni tick
    tick_pos <- seq(1, p, length.out = min(20, p))
    tick_labels <- var_names[round(tick_pos)]
    axis(1, at = tick_pos, labels = tick_labels, las = 2, cex.axis = 0.3)
    axis(2, at = tick_pos, labels = tick_labels, las = 1, cex.axis = 0.3)
  }
  
  # Aggiungi griglia per matrici piccole
  if (p <= 20) {
    abline(h = 1:p + 0.5, col = "gray", lwd = 0.5)
    abline(v = 1:p + 0.5, col = "gray", lwd = 0.5)
  }
  
  # Aggiungi scala colori
  legend_x <- par("usr")[2] + 0.02 * (par("usr")[2] - par("usr")[1])
  legend_y_bottom <- par("usr")[3]
  legend_y_top <- par("usr")[4]
  
  # Disegna barra dei colori
  legend_colors <- colorRampPalette(c("red", "white", "blue"))(100)
  legend_breaks <- seq(-1, 1, length.out = 101)
  
  for (i in 1:100) {
    rect(legend_x, legend_y_bottom + (i-1) * (legend_y_top - legend_y_bottom)/100,
         legend_x + 0.03 * (par("usr")[2] - par("usr")[1]),
         legend_y_bottom + i * (legend_y_top - legend_y_bottom)/100,
         col = legend_colors[i], border = NA)
  }
  
  # Etichette scala
  text(legend_x + 0.05 * (par("usr")[2] - par("usr")[1]), legend_y_bottom, "-1", cex = 0.7)
  text(legend_x + 0.05 * (par("usr")[2] - par("usr")[1]), legend_y_top, "1", cex = 0.7)
  text(legend_x + 0.05 * (par("usr")[2] - par("usr")[1]), 
       (legend_y_bottom + legend_y_top)/2, "0", cex = 0.7)
}
# Funzione per matrice di correlazione tra tutte le X
correlazione_tra_x <- function(data, y_col = NULL, metodo = "pearson", 
                               plot = T, titolo = "Matrice di Correlazione tra X", 
                               p_threshold = 20) {
  
  # Validazione input più flessibile
  if (!is.data.frame(data) && !inherits(data, c("tbl_df", "data.table"))) {
    stop("L'input deve essere un data.frame, tibble o data.table")
  }
  
  # Verifica esistenza y_col se specificata
  if (!is.null(y_col) && !y_col %in% names(data)) {
    warning(paste("La variabile Y specificata '", y_col, "' non esiste nel dataset"))
    y_col <- NULL
  }
  
  # Seleziona solo variabili numeriche (escludendo Y se specificata)
  var_numeriche <- sapply(data, is.numeric)
  if (!is.null(y_col)) {
    x_vars <- names(data)[var_numeriche & names(data) != y_col]
  } else {
    x_vars <- names(data)[var_numeriche]
  }
  
  if (length(x_vars) < 2) {
    stop("Servono almeno 2 variabili X numeriche per calcolare la matrice di correlazione")
  }
  
  # Seleziona dati numerici
  x_data <- data[, x_vars, drop = FALSE]
  
  # Verifica che ci siano osservazioni sufficienti
  n_complete <- sum(complete.cases(x_data))
  if (n_complete < 2) {
    stop("Non ci sono abbastanza osservazioni complete per calcolare le correlazioni")
  }
  
  # Calcola matrice di correlazione
  cor_matrix <- cor(x_data, method = metodo, use = "pairwise.complete.obs")
  
  # Verifica che la matrice sia valida
  if (any(is.na(diag(cor_matrix)))) {
    warning("Alcune variabili hanno varianza zero o problemi di calcolo")
  }
  
  # Plot se richiesto
  if (plot && !is.null(cor_matrix)) {
    p <- ncol(cor_matrix)
    
    if (p <= p_threshold) {
      # Usa corrplot per p piccolo
      plot_corrplot(cor_matrix, titolo, metodo)
    } else {
      # Usa imageplot per p grande
      plot_imageplot_matrix(cor_matrix, titolo, metodo)
    }
  }
  
  return(cor_matrix)
}

grafico_distribuzioni <- function(df) {
  # Verifica che l'input sia un dataframe, tibble o data.table
  if (!is.data.frame(df) && !inherits(df, c("tbl_df", "tbl", "data.table"))) {
    stop("L'input deve essere un dataframe, tibble o data.table")
  }
  
  # Converti in dataframe standard per compatibilità
  df <- as.data.frame(df)
  
  # Verifica che il dataframe non sia vuoto
  if (nrow(df) == 0 || ncol(df) == 0) {
    stop("Il dataframe non può essere vuoto")
  }
  
  # Carica le librerie necessarie
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork")
  }
  library(ggplot2)
  library(patchwork)
  
  plot_list <- list()
  
  for (col_name in names(df)) {
    col_data <- df[[col_name]]
    
    # Salta colonne completamente vuote
    if (all(is.na(col_data))) {
      message(paste("Saltando", col_name, ": colonna completamente vuota"))
      next
    }
    
    p <- NULL # Inizializza l'oggetto plot
    
    if (is.numeric(col_data)) {
      # Per variabili numeriche, usa istogramma
      # Rimuovi NA per il calcolo dei bin, ma mantienili nel dataset per trasparenza
      non_na_data <- col_data[!is.na(col_data)]
      
      if (length(non_na_data) == 0) {
        message(paste("Saltando", col_name, ": tutti i valori numerici sono NA"))
        next
      }
      
      # Calcola bins appropriati
      n_bins <- min(30, max(10, length(non_na_data) / 10))
      
      p <- ggplot(df, aes(x = .data[[col_name]])) +
        geom_histogram(fill = "skyblue", color = "black", bins = n_bins, alpha = 0.7) +
        labs(title = paste("Distribuzione di", col_name), 
             x = col_name, 
             y = "Frequenza",
             subtitle = if(sum(is.na(col_data)) > 0) paste("Valori NA:", sum(is.na(col_data))) else NULL) +
        theme_minimal()
      
    } else if (is.factor(col_data) || is.character(col_data) || is.logical(col_data)) {
      # Per variabili categoriche, carattere o logiche
      # Crea una copia temporanea del dataframe per evitare modifiche all'originale
      temp_df_for_plot <- df
      
      # Converti in carattere per manipolare NA
      temp_df_for_plot[[col_name]] <- as.character(temp_df_for_plot[[col_name]])
      temp_df_for_plot[[col_name]][is.na(temp_df_for_plot[[col_name]])] <- "NAs"
      
      # Ricrea il fattore includendo "NAs" come livello
      # Mantiene l'ordine originale delle categorie, mettendo NAs alla fine
      original_levels <- NULL
      if (is.factor(col_data)) {
        original_levels <- levels(col_data)
      } else {
        # Se era carattere o logico, ottieni i livelli unici non NA
        original_levels <- unique(na.omit(col_data))
      }
      
      all_levels <- c(original_levels, "NAs")
      temp_df_for_plot[[col_name]] <- factor(temp_df_for_plot[[col_name]], levels = all_levels)
      
      # Limita il numero di categorie mostrate per leggibilità
      max_categories <- 20
      if (length(levels(temp_df_for_plot[[col_name]])) > max_categories) {
        # Mantieni le top categorie (esclusi NAs) e raggruppa il resto
        freq_table <- table(temp_df_for_plot[[col_name]])
        nas_count <- freq_table["NAs"]
        non_nas_freq <- freq_table[names(freq_table) != "NAs"]
        
        if (length(non_nas_freq) > 0) {
          top_categories <- names(sort(non_nas_freq, decreasing = TRUE))[1:(max_categories-2)]
          temp_df_for_plot[[col_name]] <- as.character(temp_df_for_plot[[col_name]])
          temp_df_for_plot[[col_name]][!temp_df_for_plot[[col_name]] %in% c(top_categories, "NAs")] <- "Altri"
          
          # Riordina: top categories, Altri, NAs
          final_levels <- c(top_categories, "Altri")
          if (!is.na(nas_count) && nas_count > 0) {
            final_levels <- c(final_levels, "NAs")
          }
          temp_df_for_plot[[col_name]] <- factor(temp_df_for_plot[[col_name]], levels = final_levels)
        }
      }
      
      p <- ggplot(temp_df_for_plot, aes(x = .data[[col_name]])) +
        geom_bar(fill = "lightgreen", color = "black", alpha = 0.7) +
        labs(title = paste("Distribuzione di", col_name), 
             x = col_name, 
             y = "Conteggio") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
      
    } else {
      # Per altri tipi non supportati
      message(paste("Saltando", col_name, ": tipo di dato non supportato (", class(col_data)[1], ")"))
      next
    }
    
    if (!is.null(p)) {
      plot_list[[col_name]] <- p
    }
  }
  
  # Verifica che ci siano plot da mostrare
  if (length(plot_list) == 0) {
    message("Nessun plot è stato generato.")
    return(invisible(NULL))
  }
  
  # Disponi i plot usando patchwork
  n_plots <- length(plot_list)
  
  # Calcola layout ottimale
  if (n_plots == 1) {
    result <- plot_list[[1]]
  } else if (n_plots <= 4) {
    ncol <- 2
    result <- wrap_plots(plot_list, ncol = ncol)
  } else if (n_plots <= 9) {
    ncol <- 3
    result <- wrap_plots(plot_list, ncol = ncol)
  } else {
    # Per più di 9 plot, crea pagine multiple
    plots_per_page <- 9
    n_pages <- ceiling(n_plots / plots_per_page)
    
    message(paste("Generando", n_pages, "pagine di grafici..."))
    
    # Crea lista di pagine
    pages <- list()
    for (page in 1:n_pages) {
      start_idx <- (page - 1) * plots_per_page + 1
      end_idx <- min(page * plots_per_page, n_plots)
      page_plots <- plot_list[start_idx:end_idx]
      
      pages[[page]] <- wrap_plots(page_plots, ncol = 3)
    }
    
    # Restituisci la prima pagina e stampa le altre
    result <- pages[[1]]
    if (n_pages > 1) {
      for (i in 2:n_pages) {
        cat("\n\n")
        print(pages[[i]])
      }
    }
  }
  
  return(result)
}

plot_x_vs_y <- function(df, y_var_names,
                        numeric_x_test = c("none", "linear_model"),
                        test_type = c("parametric", "non_parametric")) {
  
  # Validazione input
  if (!is.data.frame(df) && !inherits(df, c("tbl_df", "tbl", "data.table")))
    stop("L'input deve essere un dataframe, tibble o data.table")
  df <- as.data.frame(df)
  if (nrow(df) == 0 || ncol(df) == 0) stop("Il dataframe non può essere vuoto")
  
  if (!requireNamespace("ggplot2",   quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
  if (!requireNamespace("ggpubr",    quietly = TRUE)) install.packages("ggpubr")
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
  
  numeric_x_test <- match.arg(numeric_x_test)
  test_type      <- match.arg(test_type)
  
  # Validazione variabili Y
  for (y_var_name in y_var_names) {
    if (!y_var_name %in% names(df))
      stop(paste("La variabile y_var_name (", y_var_name, ") non esiste nel dataframe."))
    if (!is.numeric(df[[y_var_name]]))
      stop(paste("La variabile y (", y_var_name, ") deve essere quantitativa (numerica)."))
  }
  
  single_y <- length(y_var_names) == 1
  x_cols   <- setdiff(names(df), y_var_names)
  
  # ── Helper: costruisce UN singolo plot (col_name ~ y_var_name) ──────────────
  make_single_plot <- function(col_name, y_var_name) {
    x_data             <- df[[col_name]]
    p                  <- NULL
    test_result_text   <- NULL
    temp_df_for_plot   <- df
    
    if (is.numeric(x_data)) {
      temp_df_clean <- na.omit(temp_df_for_plot[, c(col_name, y_var_name)])
      if (nrow(temp_df_clean) == 0) return(NULL)
      
      title_str <- if (single_y)
        paste("Relazione tra", col_name, "e", y_var_name) else y_var_name
      
      p <- ggplot(temp_df_clean, aes(x = .data[[col_name]], y = .data[[y_var_name]])) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "lm", se = TRUE, color = "blue") +
        labs(title    = title_str,
             x        = col_name,
             y        = y_var_name,
             subtitle = if (nrow(temp_df_for_plot) - nrow(temp_df_clean) > 0)
               paste("Osservazioni rimosse (NA):", nrow(temp_df_for_plot) - nrow(temp_df_clean))
             else NULL) +
        theme_minimal()
      
      if (numeric_x_test == "linear_model" && nrow(temp_df_clean) > 2) {
        tryCatch({
          lm_model  <- lm(as.formula(paste(y_var_name, "~", col_name)), data = temp_df_clean)
          model_sum <- summary(lm_model)
          if (col_name %in% rownames(model_sum$coefficients)) {
            p_value_lm <- model_sum$coefficients[col_name, "Pr(>|t|)"]
            r_squared  <- model_sum$r.squared
            test_result_text <- if (single_y)
              paste0("Modello Lineare p-value: ", format.pval(p_value_lm, digits = 3),
                     " (R² = ", round(r_squared, 3), ")")
            else
              paste0("LM p=", format.pval(p_value_lm, digits = 3),
                     " R²=", round(r_squared, 3))
          }
        }, error = function(e) {
          test_result_text <<- if (single_y) "Modello non applicabile (errore)" else "LM errore"
        })
      }
      
    } else if (is.factor(x_data) || is.character(x_data) || is.logical(x_data)) {
      temp_df_for_plot[[col_name]] <- as.character(temp_df_for_plot[[col_name]])
      temp_df_for_plot[[col_name]][is.na(temp_df_for_plot[[col_name]])] <- "NAs"
      original_levels <- if (is.factor(x_data)) levels(x_data) else unique(na.omit(x_data))
      temp_df_for_plot[[col_name]] <- factor(temp_df_for_plot[[col_name]],
                                             levels = c(original_levels, "NAs"))
      temp_df_clean <- temp_df_for_plot[!is.na(temp_df_for_plot[[y_var_name]]), ]
      if (nrow(temp_df_clean) == 0) return(NULL)
      
      # Limita categorie
      max_categories <- 15
      if (length(levels(temp_df_clean[[col_name]])) > max_categories) {
        freq_table     <- table(temp_df_clean[[col_name]])
        nas_count      <- freq_table["NAs"]
        non_nas_freq   <- freq_table[names(freq_table) != "NAs"]
        top_categories <- names(sort(non_nas_freq, decreasing = TRUE))[1:(max_categories - 2)]
        temp_df_clean[[col_name]] <- as.character(temp_df_clean[[col_name]])
        temp_df_clean[[col_name]][!temp_df_clean[[col_name]] %in% c(top_categories, "NAs")] <- "Altri"
        final_levels <- c(top_categories, "Altri")
        if (!is.na(nas_count) && nas_count > 0) final_levels <- c(final_levels, "NAs")
        temp_df_clean[[col_name]] <- factor(temp_df_clean[[col_name]], levels = final_levels)
      }
      
      title_str <- if (single_y)
        paste("Relazione tra", col_name, "e", y_var_name) else y_var_name
      
      p <- ggplot(temp_df_clean, aes(x = .data[[col_name]], y = .data[[y_var_name]])) +
        geom_boxplot(fill = "lightgreen", alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.4, color = "darkblue") +
        labs(title    = title_str,
             x        = col_name,
             y        = y_var_name,
             subtitle = if (nrow(temp_df_for_plot) - nrow(temp_df_clean) > 0)
               paste("Osservazioni rimosse (y NA):", nrow(temp_df_for_plot) - nrow(temp_df_clean))
             else NULL) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
      
      group_counts <- table(temp_df_clean[[col_name]])
      valid_groups <- names(group_counts)[group_counts >= 2]
      
      if (length(valid_groups) > 1) {
        temp_df_test <- temp_df_clean[temp_df_clean[[col_name]] %in% valid_groups, ]
        temp_df_test[[col_name]] <- droplevels(temp_df_test[[col_name]])
        num_levels <- length(levels(temp_df_test[[col_name]]))
        tryCatch({
          formula_str <- paste(y_var_name, "~", col_name)
          if (test_type == "parametric") {
            if (num_levels == 2) {
              res <- t.test(as.formula(formula_str), data = temp_df_test)
              test_result_text <- if (single_y)
                paste0("T-test p-value: ", format.pval(res$p.value, digits = 3))
              else
                paste0("t-test p=", format.pval(res$p.value, digits = 3))
            } else {
              res   <- aov(as.formula(formula_str), data = temp_df_test)
              p_val <- summary(res)[[1]]$`Pr(>F)`[1]
              test_result_text <- if (single_y)
                paste0("ANOVA p-value: ", format.pval(p_val, digits = 3))
              else
                paste0("ANOVA p=", format.pval(p_val, digits = 3))
            }
          } else {
            if (num_levels == 2) {
              res <- wilcox.test(as.formula(formula_str), data = temp_df_test)
              test_result_text <- if (single_y)
                paste0("Wilcoxon p-value: ", format.pval(res$p.value, digits = 3))
              else
                paste0("Wilcox p=", format.pval(res$p.value, digits = 3))
            } else {
              res <- kruskal.test(as.formula(formula_str), data = temp_df_test)
              test_result_text <- if (single_y)
                paste0("Kruskal-Wallis p-value: ", format.pval(res$p.value, digits = 3))
              else
                paste0("KW p=", format.pval(res$p.value, digits = 3))
            }
          }
        }, error = function(e) {
          test_result_text <<- if (single_y) "Test non applicabile (errore)" else "Test errore"
        })
      } else {
        test_result_text <- if (single_y) "Test non applicabile (<2 gruppi validi)" else "<2 gruppi validi"
      }
      
    } else {
      message(paste("Saltando", col_name, ": tipo di dato non supportato (", class(x_data)[1], ")"))
      return(NULL)
    }
    
    if (!is.null(p) && !is.null(test_result_text))
      p <- p + annotate("text", x = Inf, y = Inf,
                        label = test_result_text,
                        hjust = 1.05, vjust = 1.5,
                        size  = if (single_y) 3.5 else 3,
                        color = "red")
    p
  }
  
  # ── Modalità Y singola: comportamento originale (patchwork + return) ────────
  if (single_y) {
    y_var_name <- y_var_names[[1]]
    plot_list  <- list()
    
    for (col_name in x_cols) {
      if (length(unique(na.omit(df[[col_name]]))) <= 1) next
      p <- make_single_plot(col_name, y_var_name)
      if (!is.null(p)) plot_list[[col_name]] <- p
    }
    
    if (length(plot_list) == 0) {
      message("Nessun plot è stato generato.")
      return(invisible(NULL))
    }
    
    n_plots <- length(plot_list)
    if (n_plots == 1) {
      return(plot_list[[1]])
    } else if (n_plots <= 4) {
      return(wrap_plots(plot_list, ncol = 2))
    } else if (n_plots <= 9) {
      return(wrap_plots(plot_list, ncol = 3))
    } else {
      plots_per_page <- 9
      n_pages        <- ceiling(n_plots / plots_per_page)
      message(paste("Generando", n_pages, "pagine di grafici..."))
      pages <- lapply(seq_len(n_pages), function(page) {
        idx <- seq((page - 1) * plots_per_page + 1, min(page * plots_per_page, n_plots))
        wrap_plots(plot_list[idx], ncol = 3)
      })
      for (i in seq_along(pages)[-1]) { cat("\n\n"); print(pages[[i]]) }
      return(pages[[1]])
    }
  }
  
  # ── Modalità Y multipla: una riga per X, stampa progressiva ────────────────
  for (col_name in x_cols) {
    if (length(unique(na.omit(df[[col_name]]))) <= 1) next
    plots_for_x <- Filter(Negate(is.null),
                          lapply(y_var_names, function(y) make_single_plot(col_name, y)))
    if (length(plots_for_x) == 0) next
    
    combined <- wrap_plots(plots_for_x, nrow = 1) +
      plot_annotation(title = paste("X:", col_name),
                      theme = theme(plot.title = element_text(face = "bold", size = 13)))
    print(combined)
    cat("\n")
  }
  
  invisible(NULL)
}

plot_density_y_by_x <- function(df, y_var_names,
                                test_type    = c("parametric", "non_parametric"),
                                numeric_x_bins = c("none", "quartile")) {
  # ── Validazione input ──────────────────────────────────────────────────────
  if (!is.data.frame(df) && !inherits(df, c("tbl_df", "tbl", "data.table")))
    stop("L'input deve essere un dataframe, tibble o data.table")
  df <- as.data.frame(df)
  if (nrow(df) == 0 || ncol(df) == 0) stop("Il dataframe non può essere vuoto")
  
  for (pkg in c("ggplot2", "patchwork", "ggridges")) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  }
  library(ggplot2)
  library(patchwork)
  library(ggridges)
  
  test_type      <- match.arg(test_type)
  numeric_x_bins <- match.arg(numeric_x_bins)
  
  # Validazione variabili Y
  for (yn in y_var_names) {
    if (!yn %in% names(df))
      stop(paste("La variabile y_var_name (", yn, ") non esiste nel dataframe."))
    if (!is.numeric(df[[yn]]))
      stop(paste("La variabile y (", yn, ") deve essere numerica."))
  }
  
  single_y <- length(y_var_names) == 1
  x_cols   <- setdiff(names(df), y_var_names)
  
  # Palette discreta (fino a 15 livelli)
  PALETTE <- c(
    "#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
    "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC",
    "#D37295","#FABFD2","#8CD17D","#86BCB6","#499894"
  )
  
  # ── Helper: statistiche per caption ───────────────────────────────────────
  compute_test_text <- function(col_name, y_var_name, tmp, single) {
    grp_counts <- table(tmp[[col_name]])
    valid_grps <- names(grp_counts)[grp_counts >= 2]
    if (length(valid_grps) < 2)
      return(if (single) "Test non applicabile (<2 gruppi validi)" else "<2 gruppi")
    
    td <- tmp[tmp[[col_name]] %in% valid_grps, ]
    td[[col_name]] <- droplevels(td[[col_name]])
    nl  <- length(levels(td[[col_name]]))
    frm <- as.formula(paste(y_var_name, "~", col_name))
    
    tryCatch({
      if (test_type == "parametric") {
        if (nl == 2) {
          res <- t.test(frm, data = td)
          if (single) paste0("T-test p = ", format.pval(res$p.value, digits = 3))
          else        paste0("t=", format.pval(res$p.value, digits = 3))
        } else {
          res <- aov(frm, data = td)
          pv  <- summary(res)[[1]]$`Pr(>F)`[1]
          if (single) paste0("ANOVA p = ", format.pval(pv, digits = 3))
          else        paste0("ANOVA=", format.pval(pv, digits = 3))
        }
      } else {
        if (nl == 2) {
          res <- wilcox.test(frm, data = td)
          if (single) paste0("Wilcoxon p = ", format.pval(res$p.value, digits = 3))
          else        paste0("W=", format.pval(res$p.value, digits = 3))
        } else {
          res <- kruskal.test(frm, data = td)
          if (single) paste0("Kruskal-Wallis p = ", format.pval(res$p.value, digits = 3))
          else        paste0("KW=", format.pval(res$p.value, digits = 3))
        }
      }
    }, error = function(e) if (single) "Test non applicabile (errore)" else "errore")
  }
  
  # ── Helper: costruisce UN plot densità ────────────────────────────────────
  make_density_plot <- function(col_name, y_var_name) {
    x_raw <- df[[col_name]]
    
    # Gestione numerica: binning a quartili (opzionale)
    if (is.numeric(x_raw)) {
      if (numeric_x_bins == "none") return(NULL)
      q   <- quantile(x_raw, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
      lbl <- c("Q1","Q2","Q3","Q4")
      tmp <- df
      tmp[[col_name]] <- cut(x_raw, breaks = q, labels = lbl,
                             include.lowest = TRUE)
      if (all(is.na(tmp[[col_name]]))) return(NULL)
      bin_note <- paste0("(binned in quartili)")
    } else {
      tmp      <- df
      bin_note <- NULL
    }
    
    # Preparazione del fattore
    x_vals <- tmp[[col_name]]
    x_vals <- as.character(x_vals)
    x_vals[is.na(x_vals)] <- "NA"
    orig_levels <- if (is.factor(df[[col_name]])) levels(df[[col_name]])
    else unique(na.omit(as.character(df[[col_name]])))
    all_levels  <- c(orig_levels[orig_levels != "NA"], "NA")
    # Limita a 15 categorie
    max_cat <- 15
    if (length(all_levels) > max_cat) {
      ft   <- table(x_vals[x_vals != "NA"])
      top  <- names(sort(ft, decreasing = TRUE))[1:(max_cat - 1)]
      x_vals[!x_vals %in% c(top, "NA")] <- "Altri"
      all_levels <- c(top, "Altri")
      if (any(x_vals == "NA")) all_levels <- c(all_levels, "NA")
    }
    tmp[[col_name]] <- factor(x_vals, levels = all_levels)
    
    # Rimuovi NA in y
    tmp <- tmp[!is.na(tmp[[y_var_name]]), ]
    if (nrow(tmp) == 0) return(NULL)
    if (length(unique(tmp[[col_name]])) < 2) return(NULL)
    
    n_lvl   <- length(levels(tmp[[col_name]]))
    palette <- PALETTE[seq_len(min(n_lvl, length(PALETTE)))]
    
    # Titolo e test
    test_txt <- compute_test_text(col_name, y_var_name, tmp, single_y)
    title_str <- if (single_y)
      paste("Densità di", y_var_name, "per", col_name,
            if (!is.null(bin_note)) bin_note else "")
    else y_var_name
    
    subtitle_str <- paste0(test_txt,
                           if (!is.null(bin_note)) paste0(" — ", bin_note) else "")
    
    # Ridge plot: y = livelli di X, x = valori di Y
    # Ordine top->bottom segue l'ordine dei livelli
    tmp[[col_name]] <- factor(tmp[[col_name]],
                              levels = rev(levels(tmp[[col_name]])))
    
    p <- ggplot(tmp, aes(
      x    = .data[[y_var_name]],
      y    = .data[[col_name]],
      fill = .data[[col_name]])) +
      geom_density_ridges(alpha = 0.70, scale = 1.1,
                          color = "white", linewidth = 0.3,
                          quantile_lines = TRUE, quantiles = 2) +
      scale_fill_manual(values = palette, guide = "none") +
      labs(title    = title_str,
           subtitle = subtitle_str,
           x        = y_var_name,
           y        = col_name) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title    = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9, color = "red"),
        axis.text.y   = element_text(size = 8)
      )
    p
  }
  
  # ── Modalità Y singola ────────────────────────────────────────────────────
  if (single_y) {
    y_var_name <- y_var_names[[1]]
    plot_list  <- list()
    
    for (col_name in x_cols) {
      if (length(unique(na.omit(df[[col_name]]))) <= 1) next
      p <- make_density_plot(col_name, y_var_name)
      if (!is.null(p)) plot_list[[col_name]] <- p
    }
    
    if (length(plot_list) == 0) { message("Nessun plot generato."); return(invisible(NULL)) }
    
    n <- length(plot_list)
    if (n == 1)       return(plot_list[[1]])
    else if (n <= 4)  return(wrap_plots(plot_list, ncol = 2))
    else if (n <= 9)  return(wrap_plots(plot_list, ncol = 3))
    else {
      ppp    <- 6          # ridge plots sono alti, meglio 6 per pagina
      npages <- ceiling(n / ppp)
      message(paste("Generando", npages, "pagine..."))
      pages  <- lapply(seq_len(npages), function(pg) {
        idx <- seq((pg - 1) * ppp + 1, min(pg * ppp, n))
        wrap_plots(plot_list[idx], ncol = 2)
      })
      for (i in seq_along(pages)[-1]) { cat("\n\n"); print(pages[[i]]) }
      return(pages[[1]])
    }
  }
  
  # ── Modalità Y multipla ───────────────────────────────────────────────────
  for (col_name in x_cols) {
    if (length(unique(na.omit(df[[col_name]]))) <= 1) next
    plots_for_x <- Filter(Negate(is.null),
                          lapply(y_var_names, function(y) make_density_plot(col_name, y)))
    if (length(plots_for_x) == 0) next
    
    combined <- wrap_plots(plots_for_x, nrow = 1) +
      plot_annotation(
        title = paste("X:", col_name),
        theme = theme(plot.title = element_text(face = "bold", size = 13))
      )
    print(combined)
    cat("\n")
  }
  invisible(NULL)
}

plot_depth_profiles <- function(df,
                                y_vars,
                                color_var  = NULL,
                                bottom_col = "Bottom",
                                group_col  = "Field",
                                log_y      = FALSE,
                                ncol       = 1) {
  if (!requireNamespace("ggplot2",   quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
  library(ggplot2)
  library(patchwork)

  for (v in c(y_vars, bottom_col, group_col))
    if (!v %in% names(df)) stop(paste("Colonna non trovata:", v))
  if (!is.null(color_var) && !color_var %in% names(df))
    stop(paste("color_var non trovato:", color_var))

  col_data   <- if (!is.null(color_var)) df[[color_var]] else NULL
  is_discrete <- is.null(col_data) || is.factor(col_data) ||
                 is.character(col_data) || is.logical(col_data)
  n_levels   <- if (is_discrete && !is.null(col_data))
                  length(unique(na.omit(col_data))) else 0L

  make_panel <- function(y_var) {
    d <- df[!is.na(df[[y_var]]) & !is.na(df[[bottom_col]]), ]
    if (log_y) {
      d[[y_var]] <- log(d[[y_var]])
      y_label <- paste0("log(", y_var, ")")
    } else {
      y_label <- y_var
    }

    p <- ggplot(d, aes(x = .data[[bottom_col]], y = .data[[y_var]],
                       group = .data[[group_col]]))

    if (is.null(color_var)) {
      p <- p + geom_line(alpha = 0.5, color = "steelblue") +
               geom_point(size = 1.2, alpha = 0.7, color = "steelblue")
    } else if (is_discrete && n_levels <= 8) {
      p <- p + geom_line(aes(color = .data[[color_var]]), alpha = 0.6) +
               geom_point(aes(color = .data[[color_var]]), size = 1.2, alpha = 0.8) +
               scale_color_brewer(palette = "Set1", name = color_var)
    } else if (is_discrete) {
      p <- p + geom_line(aes(color = .data[[color_var]]), alpha = 0.5) +
               geom_point(aes(color = .data[[color_var]]), size = 1.0, alpha = 0.7) +
               scale_color_viridis_d(option = "turbo", name = color_var)
    } else {
      p <- p + geom_line(aes(color = .data[[color_var]]), alpha = 0.6) +
               geom_point(aes(color = .data[[color_var]]), size = 1.2, alpha = 0.8) +
               scale_color_viridis_c(name = color_var)
    }

    p + labs(x = paste0(bottom_col, " (cm)"), y = y_label, title = y_label) +
        theme_minimal(base_size = 11) +
        theme(legend.position = if (length(y_vars) == 1) "right" else "none")
  }

  panels <- lapply(y_vars, make_panel)

  if (length(panels) > 1 && !is.null(color_var))
    panels[[length(panels)]] <- panels[[length(panels)]] +
      theme(legend.position = "right")

  wrap_plots(panels, ncol = ncol, guides = "collect") +
    plot_annotation(
      title    = paste0("Profili verticali",
                        if (!is.null(color_var)) paste0(" — per ", color_var)),
      subtitle = if (log_y) "Scala logaritmica: linee rette ≈ decadimento esponenziale"
    )
}

#setup grafico :
colori_zone <- c(
  "Cl"     = "#B37A5C",  
  "SiCl"   = "#9F9A83",  
  "SaCl"   = "#C3A482",  
  "ClLo"   = "#9C7355",  
  "SiClLo" = "#826A4D",  
  "SaClLo" = "#B89B74",  
  "Lo"     = "#D1B993",  
  "SiLo"   = "#BCAD92",  
  "SaLo"   = "#CBAA7B",  
  "Si"     = "#A9A590",  
  "LoSa"   = "#D6C5A0",  
  "Sa"     = "#DFD4B6"   
)

# 2. Abbassa l'opacità dei colori al 50% (alpha.f = 0.5)
colori_trasparenti <- adjustcolor(colori_zone, alpha.f = 0.8)

plot_texture_triangle <- function(df,
                                   sand_col  = "PercSand",
                                   clay_col  = "PercClay",
                                   silt_col  = "PercSilt",
                                   color_var = NULL,
                                   palette   = NULL,
                                   version   = c("ggtern", "ttplot")) {
  version <- match.arg(version)

  for (col in c(sand_col, clay_col, silt_col))
    if (!col %in% names(df)) stop(paste("Colonna non trovata:", col))
  if (!is.null(color_var) && !color_var %in% names(df))
    stop(paste("color_var non trovato:", color_var))

  tex <- data.frame(SAND = df[[sand_col]], CLAY = df[[clay_col]], SILT = df[[silt_col]])

  if (!is.null(color_var)) {
    col_data <- df[[color_var]]
    lvls <- if (is.factor(col_data)) levels(col_data)
            else sort(unique(na.omit(as.character(col_data))))
    if (is.null(palette)) {
      base_pal <- c("#2166ac", "#4393c3", "#1a9850", "#74c476", "#fdae61",
                    "#d73027", "#762a83", "#e08214", "#a6761d", "#666666",
                    "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
      palette <- setNames(base_pal[seq_len(min(length(lvls), length(base_pal)))], lvls)
    }
  }

  # ── versione TT.plot ────────────────────────────────────────────────────────
  if (version == "ttplot") {
    pt_col <- if (is.null(color_var)) "black"
              else palette[as.character(df[[color_var]])]
    soiltexture::TT.plot(
      class.sys      = "USDA.TT",
      tri.data       = tex,
      main           = paste0("Soil Texture Triangle – USDA",
                              if (!is.null(color_var)) paste0("\n(", color_var, ")")),
      bg             = "white",
      frame.bg.col   = "white",
      class.p.bg.col = colori_trasparenti,
      class.line.col = "#6B4226",
      class.lab.show = "full",
      class.lab.col  = "#3A2408",
      cex.lab        = 0.85,
      pch            = 16,
      col            = pt_col,
      cex            = 1.3
    )
    if (!is.null(color_var))
      legend("topright", legend = lvls, col = palette, pch = 16,
             title = color_var, bty = "n", cex = 0.8)
    return(invisible(NULL))
  }

  # ── versione ggtern ─────────────────────────────────────────────────────────
  if (!requireNamespace("ggtern", quietly = TRUE))
    stop("Installa ggtern con: install.packages('ggtern')")
  library(ggtern)

  step    <- 0.5
  grid_bg <- expand.grid(CLAY = seq(0, 100, by = step),
                         SILT = seq(0, 100, by = step))
  grid_bg <- grid_bg[grid_bg$CLAY + grid_bg$SILT <= 100, ]
  grid_bg$SAND <- 100 - grid_bg$CLAY - grid_bg$SILT

  mat <- soiltexture::TT.points.in.classes(tri.data = grid_bg, class.sys = "USDA.TT")
  grid_bg$classe <- apply(mat, 1, function(x) {
    nm <- names(x)[x > 0]
    if (length(nm) == 0L) NA_character_ else nm[1L]
  })
  grid_bg <- grid_bg[!is.na(grid_bg$classe), ]

  centroids <- do.call(rbind, lapply(split(grid_bg, grid_bg$classe), function(d)
    data.frame(classe = d$classe[1],
               CLAY   = mean(d$CLAY),
               SILT   = mean(d$SILT),
               SAND   = mean(d$SAND))
  ))

  crop_tex <- tex
  if (!is.null(color_var)) crop_tex[[color_var]] <- df[[color_var]]

  p <- ggtern(grid_bg, aes(x = SAND, y = CLAY, z = SILT)) +
    geom_point(aes(colour = classe), shape = 15, size = 1.1, alpha = 0.85) +
    scale_colour_manual(values = colori_zone, guide = "none")

  # inherit.aes = TRUE: eredita x/y/z dal ggtern parent, evita i warning su "z"
  if (!is.null(color_var)) {
    p <- p +
      geom_point(
        data        = crop_tex,
        aes(fill    = .data[[color_var]]),
        shape       = 21, size = 2.8, stroke = 0.7, colour = "black",
        inherit.aes = TRUE
      ) +
      scale_fill_manual(values = palette, name = color_var)
  } else {
    p <- p +
      geom_point(
        data        = tex,
        shape       = 21, size = 2.8, stroke = 0.7,
        fill        = alpha("white", 0.8), colour = "black",
        inherit.aes = TRUE
      )
  }

  # Etichette: due geom_text sovrapposti (halo bianco + testo scuro)
  # position_identity() esplicito evita che ggtern rimuova il layer per PositionNudge
  p +
    geom_text(
      data        = centroids,
      aes(label   = classe),
      size        = 3.2, fontface = "bold", colour = "white",
      inherit.aes = TRUE,
      position    = position_identity()
    ) +
    geom_text(
      data        = centroids,
      aes(label   = classe),
      size        = 2.8, fontface = "bold", colour = "gray15",
      inherit.aes = TRUE,
      position    = position_identity()
    ) +
    theme_bw() +
    theme_showarrows() +
    labs(
      title = paste0("Soil Texture Triangle – USDA",
                     if (!is.null(color_var)) paste0("  (", color_var, ")")),
      x = "Sand (%)", y = "Clay (%)", z = "Silt (%)"
    ) +
    theme(
      plot.title   = element_text(face = "bold", size = 14, hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.key   = element_rect(colour = "gray80")
    )
}

