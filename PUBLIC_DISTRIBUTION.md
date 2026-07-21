# Public distribution plan for AutoRNAseq

The lowest-cost, non-coder-friendly way to publish AutoRNAseq is:

1. Keep the application as a local Docker-based tool.
2. Publish a public GitHub repository.
3. Add GitHub Releases for downloadable app bundles and startup scripts.
4. Publish this `docs/` folder through GitHub Pages as the project front door.

## User flow

1. A visitor opens the GitHub Pages site.
2. They click the latest release download link.
3. They install Docker Desktop if needed.
4. They run the startup script.
5. The app opens locally at `http://localhost:3838`.

## What to include in each release

- Source code tag
- `docker-compose.yml`
- Windows startup script
- macOS/Linux startup script
- README or quick-start PDF
- Example data if size permits

## What to document on the landing page

- What AutoRNAseq does
- That it runs locally, not as a hosted server
- The latest release link
- The one-minute startup steps
- The portfolio link
- The support email

