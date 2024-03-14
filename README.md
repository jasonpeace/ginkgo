# Backend Coding Challenge

The role of a software engineer at Ginkgo requires integration of internal and external data services as well as data persistence infrastructure. This requires versatile development skills and knowledge of secure, scalable patterns for service-to-service communication and data modeling.

In our software, we generally work with DNA and Proteins. DNA is transcribed to RNA, which then can be translated to Protein. The objective of this exercise is to create a web application that can determine whether a particular DNA strand matches any part of the DNA sequences that could encode a protein in a well-known set.

This exercise involves implementing a fairly simple sequence alignment algorithm but requires the integration of multiple components and deployment so that we can see your great work. Your solution should be closer to a minimal viable product (something to build on) than a finished production in regards to feature completeness.

It’s strongly suggested that you use the Biopython library rather than try and implement your own alignment search algorithm. Ginkgo uses Django and React (and we prefer you use them for this exercise), but if you feel more comfortable using another language or web framework, please do.

Application Overview:
- A web application serves a simple browser Javascript client that takes a DNA sequence
as a form input.
- Upon submission of a DNA sequence, an asynchronous alignment is started to find a
Protein that contains the submitted sequence.
- Proteins from the provided list of genomes are searched until one is found containing
the provided sequence. The list will be much larger in real life. Note that a single
genome may contain multiple proteins.
- The browser client allows for additional submissions before a (the first) submission has
completed.
- When a submission is complete, the web client is updated with the results, which include
the name of the protein and where the sequence was found in the protein’s sequence.
Simple text-based results is all that is required.
- A user can close their browser, and upon subsequently revisiting the web application,
their previous submissions and any active submissions are reloaded automatically.
- Users can see their previous searches and respective results.
Requirements:
- The only required technologies are Python and Javascript.
  
- A service is required for the logic; calling a command line interface via a shell invocation or invoking an external, third party web server are not desired solutions.
- Running all required processes for local development should be automated; one command. For example, if a webserver, task runner and database are required, then a single command should be able to launch all of these for the purposes of development.
We encourage questions, but please submit them in one email.
When ready for review, please submit:
1. A single Github repo containing your solution
2. The URL of the running application.

Genome List:
- NC_000852, NC_007346, NC_008724, NC_009899, NC_014637, NC_020104, NC_023423, NC_023640, NC_023719, NC_027867


## Implementation

#### The stack:
- Vite for building the React Code
- React.js for the front end implemented as a SPA
- Django for the backend
- Celery/Redis for the task runner
- Postgres for the DB

To run the code, pull the repo locally and run docker-compose build; docker-compose up.

#### Workflow
Upon loading the homepage at 127.0.0.1 or whatever IP the deployed code is runing on, the user will be taken to the home page.
The home page consists of a text area for submitting DNA queries and if there is any data, a table that shows what was entered,
its status and dates for creation/updates. There is no login for this demo. I decided implementing users would open up more
questions about user roles, what data they could see, hierarchies of permissions, etc. That can be solved later once its clear
how the tool would fit into the organization. The view methods needed to use the csrf_exempt decorator to work around this.

Once a query is submitted, a table row will appear with the data as mentioned above. This page will start polling the backend for
updates. The interval is set to 30 seconds. This is way low and would more likely be set to a much higher number in real life. I've
left it low so its easy to demonstrate. When the alignment completes, the status should flip to "COMPLETED". Before this, or after,
the user can click on the "details" button in the rightmost column. This will take them to a details page specific to this request.

On the details page, you should see the same information as the row from the home page you clicked on. If the alignment is done, there
will be addition information here; either a warning about no alignments being found or a match which shows information like protein id,
etc. If the alignment isn't done yet, a function will poll the backend until there is data one way or the other. After that it stops
polling. The interval here is the same as the home page, 30 seconds.

The alignment is done by creating a task and handing it off to celery to run. You can queue up as many requests as you like and celelry 
will work through them. When a task completes, celery hits an api endpoint on the django app to update the postgres db with results.

Biopython is used for the alignment as suggested. I've chosen to use the PairwiseAligner set in local mode. I'm not a bioinformatics
person so the implementation is probably suspect. I've set the open_gap_score to very far into the negative to reduce the number of
alignments returned. Without this, it would often find more than the len(alignments) could even return. In real life, i'd work with a
scientist or other bioinformatics people to determine the appropriate algorithm. The language in the challenge made it seem like you
wanted a match without gaps or mismatches so that is what i've attempted to implement. Its probably not the best use of the aligner.
The code for this can be found in app/protein_app/modules/alignment.py. The genomes specified in the challenge were downloaded 
in a genbank format and stored in their own folder under /app. When the aligner runs, it uses SeqIO to parses these and looks for 
features that are proteins. Their dna is extracted from the origin/parent and used as the target for the alignment.

#### Vite
Vite is being used to compile React components into the static dir used by Django. This is configured for more of a SPA style 
deployment where Django will host the page that brings up the SPA and can otherwise serve api endpoints. Note that the typescript
precompiler is not in strict mode. That can be changed in the tsconfig.json found in the vite-project dir.

You can run "npm run dev". This will start the vite dev server up in watch mode.

There are two other notable commands: build_dev and build. 

The 'build_dev' command will process the files via typescript and bundle them in the dev static dir where django can find them.
Sourcemaps are created inline. The 'build' command also runs the typescript preprocessor but then minifies them and bundles them 
for production.

#### Improvements

As the challenge mentions, this is a place to start from which we could improve. There are obvious next steps.

- Add users and permissions 
- Tune the alignment algorithm to be more useful
- expand the interface to allow configuring the aligner properties while also adding some guardrails to keep the alignment numbers
reasonable
- Improve the dev environment by using external volumes for the code rather than copying via the build. This would allow django to 
see the python changes but to also allow vite to run in watch mode so react changes could be served instantly
- The UI should match other apps used at the company










