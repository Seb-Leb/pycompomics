import os
from os.path import join as opj
import subprocess
import yaml
from shutil import copyfile

report_ids = {
    'Certificate of Analysis': '0',
    'Default Hierarchical Report': '1',
    'Default PSM Phosphorylation Report': '2',
    'Default PSM Report': '3',
    'Default PSM Report with non-validated matches': '4',
    'Default Peptide Phosphorylation Report': '5',
    'Default Peptide Report': '6',
    'Default Peptide Report with non-validated matches': '7',
    'Default Protein Phosphorylation Report': '8',
    'Default Protein Report': '9',
    'Default Protein Report with non-validated matches': '10',
    'Extended PSM Report': '11',
    'Extended PSM Annotation Report': '12'
             }

# TODO: implement logging

class SearchGUI:
    def __init__(self, fasta_db, mgf_path, out_dir, exp_name, compomics_path, searchgui_version, db_cache, ptm_config_json,
            n_threads=4, tmp_dir='tmp/', protein_fdr=1., ms_level='high', generate_fasta_db=False):
        self.fasta_db         = fasta_db    
        self.mgf_path         = mgf_path
        self.out_dir          = out_dir
        self.exp_name         = exp_name
        self.compomics_path   = compomics_path
        self.searchgui_path   = opj(compomics_path, searchgui_version, f'{searchgui_version}.jar')
        self.parameters_cache = opj(compomics_path, 'params_cache')
        self.db_cache         = db_cache
        self.ptm_config       = opj(compomics_path, 'ptm', ptm_config_json)
        self.tmp_dir          = tmp_dir
        self.n_threads        = n_threads
        if not os.path.exists(self.tmp_dir):
            os.mkdir(self.tmp_dir)
        
        db_path = fasta_db.split('.')[0] + '_concatenated_target_decoy.fasta'
        if os.path.isfile(db_path):
            self.db_path = db_path
        elif generate_fasta_db:
            self.gen_fastadb_decoy()

        self.report_dir   = opj(self.out_dir, 'reports')
        self.log_dir      = opj(self.out_dir, 'logs')
        self.protein_fdr  = protein_fdr
        self.ms_level     = ms_level
        
        if not os.path.exists(self.report_dir):
            os.mkdir(self.report_dir)
        

    def run_search(self, ):
        search_cmd = self.get_search_cmd()
        result = subprocess.run(search_cmd.split(), capture_output=True)
        print(result.stdout.decode())
        print(result.stderr.decode())

    def set_search_params(self, **kwargs):
        with open('./pycompomics/searchgui_default_params.yml', 'r') as f:
            def_params = yaml.full_load(f)
        params = def_params['ms-common']
        for p in params:
            if type(params[p]) == list:
                params[p] = '"{}"'.format(',&'.join(params[p]))

        for k in kwargs:
            params[k] = kwargs[k]
        params.update(def_params['ms-level'][self.ms_level])
        self.search_engines = def_params['search-engines']

        self.id_params_path = opj(self.out_dir, 'id_params.par')
        params['out'] = self.id_params_path
        params['db']  = self.db_path
        self.params = params

        cmd = f'java -cp {self.searchgui_path} eu.isas.searchgui.cmd.IdentificationParametersCLI '
        for k, v in self.params.items():
            cmd += f'-{k} {v} '
        cmd += f'-ptm_configuration {self.ptm_config} '
        cmd = [c.replace('&', ' ') for c in cmd.split()]
        result = subprocess.run(cmd, capture_output=True)
        print(result.stdout.decode())
        print(result.stderr.decode())

    def gen_fastadb_decoy(self):
        print('Generating DB with decoys.')
        self.db_path = opj(self.db_cache, self.fasta_db.replace('.fasta', '_concatenated_target_decoy.fasta'))
        # TODO load cache...
        cmd = f'java -cp {self.searchgui_path} eu.isas.searchgui.cmd.FastaCLI -in {self.fasta_db} -decoy'
        result = subprocess.run(cmd.split(), capture_output=True)
        print(result.stdout.decode())
        print(result.stderr.decode())

    def get_search_cmd(self, **kwargs):
        cmd  = f'java -Xmx27G -cp {self.searchgui_path} eu.isas.searchgui.cmd.SearchCLI '
        cmd += f'-spectrum_files {self.mgf_path} '
        cmd += f'-output_folder {self.out_dir}/search_results/ '
        cmd += f'-id_params {self.id_params_path} '
        for k,v in self.search_engines.items():
            cmd += f'-{k} {v} '
        cmd += f'-protein_fdr {self.protein_fdr} '
        cmd += f'-threads {self.n_threads} '
        return cmd

class PeptideShaker:
    def __init__(self, searchgui, peptideshaker_version, fasta, compomics_path=None, sample_name = 'sample', replicate = 1):
        self.searchgui = searchgui
        self.tmp_dir = searchgui.tmp_dir
        self.out_dir = searchgui.out_dir
        self.compomics_path = compomics_path if compomics_path is not None else sg.compomics_path
        self.peptideshaker_path = opj(self.compomics_path, peptideshaker_version, f'{peptideshaker_version}.jar')
        self.sample_name = sample_name
        self.replicate   = replicate
        self.fasta = fasta

        cmd = f'java -cp {self.peptideshaker_path} eu.isas.peptideshaker.cmd.PathSettingsCLI -temp_folder {self.tmp_dir} -ptm_configuration {self.searchgui.ptm_config}'
        result = subprocess.run(cmd, capture_output=True, shell=True)
        print(result.stdout.decode())
        print(result.stderr.decode())

    def run(self, ):
        cmd = f'java -Xmx27G -cp {self.peptideshaker_path} eu.isas.peptideshaker.cmd.PeptideShakerCLI '
        cmd += f'-temp_folder {self.tmp_dir} '
        cmd += f'-experiment {self.searchgui.exp_name} '
        cmd += f'-sample {self.sample_name} '
        cmd += f'-replicate {self.replicate} '
        cmd += f'-fasta_file {self.fasta} '
        #cmd += f'-id_params {self.searchgui.id_params_path} '
        cmd += f'-identification_files {self.searchgui.out_dir}/search_results/searchgui_out.zip '
        cmd += f'-spectrum_files {self.searchgui.mgf_path} '
        cmd += f'-out {self.out_dir}/{self.searchgui.exp_name}.cpsx '
        cmd += f'-zip {self.out_dir}/{self.searchgui.exp_name}.cpsx.zip '
        result = subprocess.run(cmd, capture_output=True, shell=True)
        print(result.stdout.decode())
        print(result.stderr.decode())

    def generate_reports(self, reports=[str(x) for x in range(9)]):
        if type(reports) is not list:
            raise('reports argument must be of type list.')

        self.out_reports_dir = opj(self.out_dir, 'peptideshaker_reports')
        if not os.path.exists(self.out_reports_dir):
            os.mkdir(self.out_reports_dir)

        unavailable_reports = [x for x in reports if not str(x).isdigit() and x not in report_ids]
        if unavailable_reports:
            print(f'these reports are not available: {unavailable_reports}')
        reports = [x if str(x).isdigit() else report_ids[x] for x in reports]
        reports_n = ', '.join([str(x) for x in reports])

        cmd  = f'java -cp {self.peptideshaker_path} eu.isas.peptideshaker.cmd.ReportCLI '
        cmd += f'-in {self.out_dir}/{self.searchgui.exp_name}.cpsx.zip '
        cmd += f'-out_reports {self.out_reports_dir} '
        cmd += f'-reports "{reports_n}"'
        result = subprocess.run(cmd, capture_output=True, shell=True)
        print(result.stdout.decode())
        print(result.stderr.decode())
