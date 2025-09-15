import re
import json
from pathlib import Path
from collections import defaultdict

class ProcessDependencyAnalyzer:
    """Fixed analyzer with robust Nextflow parsing"""
    
    def __init__(self, workflow_dir='.', schema_file=None, debug=True):
        self.workflow_dir = Path(workflow_dir)
        self.debug = debug
        self.processes = {}
        self.workflows = {}
        self.schema_fields = {}
        
        # Load schema
        if schema_file:
            self.load_schema(schema_file)
        else:
            self.auto_load_schema()
    
    def auto_load_schema(self):
        """Auto-find schema file"""
        for path in [
            self.workflow_dir / 'gos-assets' / 'nf-gos' / 'assets' / 'schema_input.json',
            self.workflow_dir / 'gos-assets' / 'nf-gos' / 'assets' / 'schema_inputs.json',
            self.workflow_dir / 'assets' / 'schema_input.json',
            self.workflow_dir / 'assets' / 'schema_inputs.json'
        ]:
            if path.exists():
                self.load_schema(path)
                break
    
    def load_schema(self, schema_file):
        """Load schema fields"""
        try:
            with open(schema_file, 'r') as f:
                schema = json.load(f)
            
            if 'items' in schema and 'properties' in schema['items']:
                for field_name in schema['items']['properties']:
                    self.schema_fields[field_name] = True
            
            if self.debug:
                print(f"Loaded {len(self.schema_fields)} schema fields")
        except Exception as e:
            print(f"Warning: Could not load schema: {e}")
    
    def analyze_all(self):
        """Main analysis entry point"""
        for nf_file in self.workflow_dir.rglob('*.nf'):
            if 'test' not in str(nf_file).lower() and 'work/' not in str(nf_file):
                self.extract_processes(nf_file)
                self.extract_workflows(nf_file)
    
    def extract_processes(self, file_path):
        """Extract processes from file with better brace handling"""
        try:
            content = file_path.read_text()
        except:
            return
        
        # Find process starts
        process_starts = []
        for match in re.finditer(r'process\s+(\w+)\s*\{', content):
            process_starts.append((match.start(), match.group(1), match.end()))
        
        # For each process, find its closing brace
        for i, (start, name, body_start) in enumerate(process_starts):
            # Find matching closing brace
            brace_count = 1
            pos = body_start
            
            while pos < len(content) and brace_count > 0:
                if content[pos] == '{':
                    brace_count += 1
                elif content[pos] == '}':
                    brace_count -= 1
                pos += 1
            
            if brace_count == 0:
                process_body = content[body_start:pos-1]
                
                if self.debug:
                    print(f"\nAnalyzing process: {name}")
                    # Show first part of body
                    preview = process_body[:300].replace('\n', '\\n')
                    print(f"  Body preview: {preview}...")
                
                # Extract components
                inputs = self.extract_inputs(process_body, name)
                outputs = self.extract_outputs(process_body, name)
                
                self.processes[name] = {
                    'file': str(file_path.relative_to(self.workflow_dir)),
                    'inputs': inputs,
                    'outputs': outputs,
                    'conditions': self.extract_when(process_body),
                    'directives': self.extract_directives(process_body),
                    'dependencies': {
                        'requires_fields': self.map_inputs_to_schema(inputs, name),
                        'produces_fields': self.map_outputs_to_schema(outputs, name),
                        'depends_on': []
                    }
                }
    
    def extract_inputs(self, body, process_name):
        """Extract process inputs with better parsing"""
        result = {'tuples': [], 'files': [], 'values': []}
        
        # Find input block - handle both spaces and tabs
        input_match = re.search(r'input:[\s\t]*(.*?)(?=[\s\t]*output:|[\s\t]*script:|[\s\t]*exec:|[\s\t]*when:|$)', 
                               body, re.DOTALL | re.IGNORECASE)
        if not input_match:
            if self.debug:
                print(f"  No input block found for {process_name}")
            return result
        
        input_block = input_match.group(1)
        if self.debug:
            print(f"  Input block: {input_block[:100]}...")
        
        # Handle tuple inputs - match the entire tuple declaration
        # This handles multi-line tuples
        tuple_pattern = r'tuple\s+(.*?)(?=\n[\s\t]*(?:tuple|path|file|val|output|script|exec|when)|$)'
        for match in re.finditer(tuple_pattern, input_block, re.DOTALL):
            tuple_content = match.group(1).strip()
            if self.debug:
                print(f"  Found tuple: {tuple_content}")
            
            components = []
            # Extract val(xxx), path(xxx), file(xxx) from tuple
            for comp_match in re.finditer(r'(val|path|file)\s*\(([^)]+)\)', tuple_content):
                comp_type = comp_match.group(1)
                comp_name = comp_match.group(2).strip().strip('"').strip("'")
                components.append({'type': comp_type, 'name': comp_name})
                
                # Track files
                if comp_type in ['path', 'file']:
                    result['files'].append({
                        'name': comp_name,
                        'type': comp_type,
                        'from_tuple': True
                    })
            
            if components:
                result['tuples'].append(components)
        
        # Handle standalone path/file inputs
        for match in re.finditer(r'^[\s\t]*(path|file)\s*\(([^)]+)\)', input_block, re.MULTILINE):
            if 'tuple' not in match.string[max(0, match.start()-20):match.start()]:
                result['files'].append({
                    'name': match.group(2).strip().strip('"').strip("'"),
                    'type': match.group(1),
                    'from_tuple': False
                })
        
        # Handle val inputs
        for match in re.finditer(r'^[\s\t]*val\s*\(([^)]+)\)', input_block, re.MULTILINE):
            if 'tuple' not in match.string[max(0, match.start()-20):match.start()]:
                result['values'].append(match.group(1).strip())
        
        return result
    
    def extract_outputs(self, body, process_name):
        """Extract process outputs with better parsing"""
        result = {'tuples': [], 'files': [], 'emits': {}}
        
        # Find output block - handle both spaces and tabs
        output_match = re.search(r'output:[\s\t]*(.*?)(?=[\s\t]*script:|[\s\t]*exec:|[\s\t]*when:|$)', 
                                body, re.DOTALL | re.IGNORECASE)
        if not output_match:
            if self.debug:
                print(f"  No output block found for {process_name}")
            return result
        
        output_block = output_match.group(1)
        if self.debug:
            print(f"  Output block: {output_block[:100]}...")
        
        # Handle outputs with emit - capture the whole declaration
        emit_pattern = r'(tuple|path|file)\s+(.*?)\s*,?\s*emit:\s*(\w+)'
        for match in re.finditer(emit_pattern, output_block, re.DOTALL):
            output_type = match.group(1)
            content = match.group(2).strip()
            emit_name = match.group(3)
            
            result['emits'][emit_name] = content
            
            if output_type == 'tuple':
                components = []
                for comp_match in re.finditer(r'(val|path|file)\s*\(([^)]+)\)', content):
                    comp_type = comp_match.group(1)
                    comp_name = comp_match.group(2).strip().strip('"').strip("'")
                    components.append({'type': comp_type, 'name': comp_name})
                    
                    if comp_type in ['path', 'file']:
                        result['files'].append({
                            'pattern': comp_name,
                            'type': comp_type,
                            'emit': emit_name
                        })
                
                if components:
                    result['tuples'].append(components)
            else:
                # Single path/file output
                pattern_match = re.search(r'\(([^)]+)\)', content)
                if pattern_match:
                    result['files'].append({
                        'pattern': pattern_match.group(1).strip().strip('"').strip("'"),
                        'type': output_type,
                        'emit': emit_name
                    })
        
        return result
    
    def map_inputs_to_schema(self, inputs, process_name):
        """Map process inputs to schema fields based on actual input patterns"""
        fields = []
        
        # Check files from inputs
        for file_input in inputs.get('files', []):
            name = file_input['name']
            
            # Check each schema field to see if this input matches it
            for schema_field, _ in self.schema_fields.items():
                schema_lower = schema_field.lower()
                name_lower = name.lower()
                
                # Direct pattern matching
                if 'vcf' in name_lower and 'vcf' in schema_lower:
                    # More specific VCF matching
                    if 'raw' in name_lower and 'raw' in schema_lower:
                        fields.append(schema_field)
                    elif 'raw' not in name_lower and schema_field == 'vcf':
                        fields.append(schema_field)
                elif 'bam' in name_lower and schema_field == 'bam':
                    fields.append('bam')
                elif 'cram' in name_lower and schema_field == 'cram':
                    fields.append('cram')
                elif 'fastq' in name_lower or 'fq' in name_lower or 'reads' in name_lower:
                    if schema_field in ['fastq_1', 'fastq_2']:
                        fields.append(schema_field)
                elif 'hets' in name_lower and schema_field == 'hets':
                    fields.append('hets')
                elif 'cov' in name_lower and 'cov' in schema_lower:
                    fields.append(schema_field)
                elif 'seg' in name_lower and 'seg' in schema_lower:
                    fields.append(schema_field)
                elif name == schema_field:
                    # Direct name match
                    fields.append(schema_field)
        
        return list(set(fields))
    
    def map_outputs_to_schema(self, outputs, process_name):
        """Map process outputs to schema fields"""
        fields = []
        
        # Check file outputs
        for file_output in outputs.get('files', []):
            pattern = file_output.get('pattern', '').lower()
            
            # Map based on output patterns
            if 'ffpe_filtered' in pattern or 'chimera' in process_name.lower():
                if 'RAW' in process_name:
                    if 'structural_variants_raw_chimera_filtered' in self.schema_fields:
                        fields.append('structural_variants_raw_chimera_filtered')
                else:
                    if 'structural_variants_chimera_filtered' in self.schema_fields:
                        fields.append('structural_variants_chimera_filtered')
            elif '.bam' in pattern or '*.bam' in pattern or pattern.endswith('.bam'):
                # Handle BAM outputs - including ${prefix}.bam patterns
                if 'bam' in self.schema_fields:
                    fields.append('bam')
            elif '.cram' in pattern or '*.cram' in pattern or pattern.endswith('.cram'):
                if 'cram' in self.schema_fields:
                    fields.append('cram')
            elif 'vcf' in pattern:
                # Generic VCF output
                if 'raw' in pattern and 'vcf_raw' in self.schema_fields:
                    fields.append('vcf_raw')
                elif 'vcf' in self.schema_fields:
                    fields.append('vcf')
        
        # Check process name for hints
        if 'SAMTOOLS_MERGE' in process_name and 'bam' in self.schema_fields:
            # SAMTOOLS_MERGE always produces BAM files
            if 'bam' not in fields:
                fields.append('bam')
        
        return list(set(fields))
    
    def extract_when(self, body):
        """Extract when conditions"""
        match = re.search(r'when:[\s\t]*\n(.*?)(?=\n[\s\t]*(?:input|output|script|exec):|$)', 
                         body, re.DOTALL)
        return match.group(1).strip() if match else None
    
    def extract_directives(self, body):
        """Extract process directives"""
        directives = {}
        
        # Handle multi-line container directive
        container_pattern = r'container\s+"?\$\{[^}]+\}\s*\?\s*[\'"]([^\'":]+)[\'"]?\s*:\s*[\'"]([^\'":]+)[\'"]?\s*\}'
        container_match = re.search(container_pattern, body, re.DOTALL)
        if container_match:
            # Extract both the singularity and docker container names
            sing_container = container_match.group(1).replace('docker://', '')
            docker_container = container_match.group(2)
            # Use the docker container (they're usually the same)
            directives['container'] = docker_container
        else:
            # Fallback to simple pattern
            simple_container = re.search(r'container\s+[\'"]([^\'"\n]+)[\'"]', body)
            if simple_container:
                directives['container'] = simple_container.group(1)
        
        # Other directives
        patterns = {
            'label': r'label\s+[\'"]?([^\'"]+)[\'"]?',
            'tag': r'tag\s+[\'"]?([^\'"]+)[\'"]?',
            'cpus': r'cpus\s+(\d+)',
            'memory': r'memory\s+[\'"]?([^\'"]+)[\'"]?',
            'time': r'time\s+[\'"]?([^\'"]+)[\'"]?'
        }
        
        for key, pattern in patterns.items():
            match = re.search(pattern, body)
            if match:
                directives[key] = match.group(1)
        
        return directives
    
    def extract_workflows(self, file_path):
        """Extract workflows"""
        try:
            content = file_path.read_text()
        except:
            return
        
        workflow_pattern = r'workflow\s+(\w+)\s*\{([^}]+(?:\{[^}]*\}[^}]*)*)\}'
        
        for match in re.finditer(workflow_pattern, content, re.DOTALL):
            workflow_name = match.group(1)
            workflow_body = match.group(2)
            
            self.workflows[workflow_name] = {
                'file': str(file_path.relative_to(self.workflow_dir)),
                'calls': self.extract_workflow_calls(workflow_body)
            }
    
    def extract_workflow_calls(self, body):
        """Extract process calls from workflow"""
        calls = []
        pattern = r'\b([A-Z][A-Z0-9_]*)\s*\('
        for match in re.finditer(pattern, body):
            name = match.group(1)
            if name not in ['Channel', 'EMPTY', 'VALUE', 'PATH']:
                calls.append(name)
        return calls
    
    def generate_dependency_config(self):
        """Generate output config"""
        return {
            'processes': {
                name: {
                    'file': info['file'],
                    'requires': info['dependencies']['requires_fields'],
                    'produces': info['dependencies']['produces_fields'],
                    'depends_on': info['dependencies']['depends_on'],
                    'conditions': info['conditions'],
                    'resources': info['directives']
                }
                for name, info in self.processes.items()
            },
            'workflows': self.workflows,
            'schema_fields': list(self.schema_fields.keys())
        }
    
    def print_summary(self, show_tree=True, show_groups=False):
        """Print summary"""
        print(f"\n=== Pipeline Analysis ===")
        print(f"Processes: {len(self.processes)}")
        print(f"Workflows: {len(self.workflows)}")
        print(f"Schema fields: {len(self.schema_fields)}")
        
        # Debug specific process
        if 'SV_CHIMERA_FILTER' in self.processes:
            print("\n=== SV_CHIMERA_FILTER Debug ===")
            proc = self.processes['SV_CHIMERA_FILTER']
            print(f"Inputs: {proc['inputs']}")
            print(f"Outputs: {proc['outputs']}")
            print(f"Required fields: {proc['dependencies']['requires_fields']}")
            print(f"Produced fields: {proc['dependencies']['produces_fields']}")