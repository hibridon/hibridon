pipeline {
    agent {label 'physix_agent'}
    environment {
        HIBRIDON_ROOT_PATH = '${PWD}'
    }
    stages {
        stage('Building hibridon...') {
            steps {
                sh 'make HIBRIDON_ROOT_PATH=$HIBRIDON_ROOT_PATH build'
            }
        }
        stage('Testing hibridon...') {
            steps {
                sh 'make HIBRIDON_ROOT_PATH=$HIBRIDON_ROOT_PATH test'
            }
        }
        stage('Cleaning up...') {
            steps {
                echo 'Cleaning workspace...'
                sh 'git clean -d -x -f'
                // will remove untracked files, including directories (-d) and files ignored by git (-x). Replace the -f argument with -n to perform a dry-run or -i for interactive mode, and it will tell you what will be removed.'
            }
        }        
    }
    post
    {
        // always, success, failure, unstable, changed
        success
        {
            // mail bcc: '', body: "<b>Example</b><br>Project: ${env.JOB_NAME} <br>Build Number: ${env.BUILD_NUMBER} <br> URL de build: ${env.BUILD_URL}", cc: '', charset: 'UTF-8', from: '', mimeType: 'text/html', replyTo: '', subject: "ERROR CI: Project name -> ${env.JOB_NAME}", to: "foo@foomail.com";
            mail bcc: '', body: "<b>Example</b><br>Project: ${env.JOB_NAME} <br>Build Number: ${env.BUILD_NUMBER} <br>Build URL: ${env.BUILD_URL}", cc: '', charset: 'UTF-8', from: '', mimeType: 'text/html', replyTo: '', subject: "CI build succeeded for ${env.JOB_NAME}", to: "guillaume.raffy@univ-rennes1.fr";
        }
        failure
        {
            // mail bcc: '', body: "<b>Example</b><br>Project: ${env.JOB_NAME} <br>Build Number: ${env.BUILD_NUMBER} <br> URL de build: ${env.BUILD_URL}", cc: '', charset: 'UTF-8', from: '', mimeType: 'text/html', replyTo: '', subject: "ERROR CI: Project name -> ${env.JOB_NAME}", to: "foo@foomail.com";
            mail bcc: '', body: "<b>Example</b><br>Project: ${env.JOB_NAME} <br>Build Number: ${env.BUILD_NUMBER} <br>Build URL: ${env.BUILD_URL}", cc: '', charset: 'UTF-8', from: '', mimeType: 'text/html', replyTo: '', subject: "CI build failed for ${env.JOB_NAME}", to: "guillaume.raffy@univ-rennes1.fr, benjamin.desrousseaux@univ-rennes1.fr, francois.lique@univ-rennes1.fr";
        }
    }
}
