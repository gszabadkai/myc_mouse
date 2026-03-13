# rebuild_all.R

message("🔁 Starting pipeline rebuild...")

# List all scripts in order
scripts <- list.files("scripts", pattern = "^\\d{2}_.*\\.R$", full.names = TRUE)
scripts <- sort(scripts)

# Run each script
for (script in scripts) {
  message("\n📄 Running: ", script)
  tryCatch(
    source(script),
    error = function(e) {
      message("❌ Error in ", script, ":\n", e$message)
      stop("Pipeline aborted.")
    }
  )
}

message("\n✅ Pipeline completed successfully.")
