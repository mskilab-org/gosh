import shutil

def get_environment_defaults():
    import socket
    nyu_defaults = {
        'pipeline-dir': "/gpfs/data/imielinskilab/git/mskilab/nf-gos",
        'profile': "nyu",
        'nextflow_module': "nextflow/23.04.4",
        'aws_module': 'aws-cli',
        'JAVA_HOME': '/gpfs/share/apps/jdk/17u028',
        'datasets_json': '/gpfs/data/imielinskilab/external/mskiweb/xanthc01/case-report/datasets.json'
    }
    nygc_defaults = {
        'pipeline-dir': "/gpfs/commons/groups/imielinski_lab/git/nf-gos",
        'profile': "nygc",
        'nextflow_module': "nextflow/23.10.0",
        'aws_module': 'awscli',
        'JAVA_HOME': '/nfs/sw/java/jdk-17.0.4'
    }
    mapping = {
        'fn-': nyu_defaults,
        'a100-': nyu_defaults,
        'a40-': nyu_defaults,
        'cn-': nyu_defaults,
        'gn-': nyu_defaults,
        'gpu-': nyu_defaults,
        'bigpurple': nyu_defaults,

        'mskilab0': nygc_defaults,
        'pe2': nygc_defaults,
        'ne1': nygc_defaults,
    }
    hostname = socket.gethostname()
    for prefix, defaults in mapping.items():
        if hostname.startswith(prefix):
            print(f"Detected environment: {hostname}")
            print(f"Using defaults for {hostname}")
            return defaults
    return {}

def load_required_modules(env_defaults):
    import os
    required_commands = ['nextflow', 'aws', 'singularity']
    modules_to_load = []
    load_modules_command = ""

    # Check for 'nextflow' command
    if shutil.which('nextflow') is None:
        nextflow_module = env_defaults.get('nextflow_module', 'nextflow')
        if not nextflow_module:
            nextflow_module = 'nextflow'
        modules_to_load.append(nextflow_module)
        print(f"'nextflow' command not found. Loading module '{nextflow_module}'.")
    else:
        print("'nextflow' command is already available.")

    # Check for 'aws' command
    if shutil.which('aws') is None:
        modules_to_load.append(env_defaults.get('aws_module', 'aws-cli'))
        print(f"'aws' command not found. Loading module '{env_defaults.get('aws_module', 'aws-cli')}'.")
    else:
        print("'aws' command is already available.")

    # Check for 'singularity' command
    if shutil.which('singularity') is None:
        modules_to_load.append('singularity')
        print("'singularity' command not found. Loading module 'singularity'.")
    else:
        print("'singularity' command is already available.")

    # Check if JAVA_HOME is set to correct path
    default_java = env_defaults.get('JAVA_HOME')
    envvar_java = os.environ.get("JAVA_HOME")
    if 'JAVA_HOME' in os.environ and not default_java:
        print(f"JAVA_HOME is set to {envvar_java}")
    else:
        if default_java != envvar_java:
            print(f"Setting JAVA_HOME to {default_java}")
            os.environ['JAVA_HOME'] = default_java
        else:
            print(f"JAVA_HOME is already set to {default_java}")

    # Build the load modules command string
    for module in modules_to_load:
        load_modules_command += f"module load {module} && "

    return load_modules_command
