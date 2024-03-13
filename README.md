# ginkgo


TODOs:
add typescript hinting to react code
deploy to aws the easy way
deploy to aws using github

hot loading of js would be nice
hot loading of python would be nice








# vite
Vite is being used to compile react components into the static dir used by Django. This is configued for more of a SPA style deployment where django will host the page that brings up the SPA and can otherwise serve api endpoints. Note that the typescript precompiler is not in strict mode. That can be changed in the tsconfig.json found in the vite-project dir.

vite was installed by running:

		npm create vite@latest in the /app directory
			use vite-project as the dir name.
			use react
			use typescript
		Change into vite-project
		Run npm install


Run npm run dev. This will start the vite dev server up in watch mode.

There are two other notable commands: build_dev and build. 

The 'build_dev' command will process the files via typescript and bundle them in the dev static dir where django can find it.
The 'build' command also runs the typescript preprocessor but then minifies them and bundles them for production.
