# ginkgo


TODOs:
add typescript hinting to react code
look at material ui
repush code of course

deploy to aws using github

hot loading of js would be nice
hot loading of python would be nice

what happens if the celery task doesn't find a match? Test this case.

# Postgres

I'm using Postgres as the database. Its configured to ahve a persitant volume so data won't disappear each time the container is torn down. Note that sqlite is being used for the tests.

# Users

For this demo, I did not implement login and user accounts. The view methods needed to use the csrf_exempt decorator to work around this. 


# Django Containerization

Django currently has its code copied into the container during the docker-compose build. An improvment would be to alter the dockerfile and compose to use an external mount so it can see changes. This would certainly speed up React development.


# Celery/Redis 

Celery is used in a very simple manner to allow for processes to run outside the scope of the http request. The task calls an api to send back its data.  


# Vite
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

The 'build_dev' command will process the files via typescript and bundle them in the dev static dir where django can find them. Sourcemaps are created inline.
The 'build' command also runs the typescript preprocessor but then minifies them and bundles them for production.
