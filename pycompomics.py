import os
from os.path import join as opj
import subprocess
import yaml
from shutil import copyfile

peptideshaker_path  = './tools/compomics/PeptideShaker-1.16.42/PeptideShaker-1.16.42.jar'

class SearchGUI:
    def __init__(self, fasta_db, mgf_path, out_dir, exp_name, compomics_path, searchgui_version, db_cache, ptm_config_json,
            tmp_dir='tmp/', protein_fdr=1., ms_level='high'):
        self.fasta_db         = fasta_db    
        self.mgf_path         = mgf_path
        self.out_dir          = out_dir
        self.exp_name         = exp_name
        self.searchgui_path   = opj(compomics_path, searchgui_version, f'{searchgui_version}.jar')
        self.parameters_cache = opj(compomics_path, 'params_cache')
        self.db_cache         = db_cache
        self.ptm_config       = opj(compomics_path, 'ptm', ptm_config_json)
        self.tmp_dir          = tmp_dir
        if not os.path.exists(self.tmp_dir):
            os.mkdir(self.tmp_dir)
        
        db_path = fasta_db.split('.')[0] + '_concatenated_target_decoy.fasta'
        if os.path.isfile(db_path):
            self.db_path = db_path
        else:
            self.gen_fastadb_decoy()

        self.report_dir   = os.path.join(self.out_dir, 'reports')
        self.log_dir      = os.path.join(self.out_dir, 'logs')
        self.ms_level     = ms_level
        self.protein_fdr  = protein_fdr
        if not os.path.exists(self.report_dir):
            os.mkdir(self.report_dir)

    def run_search(self, ):
        search_cmd = self.get_search_cmd()
        result = subprocess.run(search_cmd.split(), capture_output=True)
        print(result.stdout.decode())
        print(result.stderr.decode())

    def set_search_params(self, **kwargs):
        with open('./pycompomics/pycompomics/searchgui_default_params.yml', 'r') as f:
            def_params = yaml.full_load(f)
        params = def_params['ms-common']
        for p in params:
            if type(params[p]) == list:
                params[p] = '"{}"'.format(',&'.join(params[p]))

        for k in kwargs:
            params[k] = kwargs[k]
        params.update(def_params['ms-level'][self.ms_level])
        self.search_engines = def_params['search-engines']

        self.id_params_path = os.path.join(self.out_dir, 'id_params.par')
        params['out'] = self.id_params_path
        params['db']  = self.db_path
        self.params = params

        cmd = 'java -cp {} eu.isas.searchgui.cmd.IdentificationParametersCLI '.format(self.searchgui_path)
        for k, v in self.params.items():
            cmd += '-{} {} '.format(k, v)
        cmd += '-ptm_configuration {} '.format(self.ptm_config)
        cmd = [c.replace('&', ' ') for c in cmd.split()]
        result = subprocess.run(cmd, capture_output=True)
        print(result.stdout.decode())
        print(result.stderr.decode())

    def gen_fastadb_decoy(self):
        print('Generating DB with decoys.')
        cmd = 'java -cp {} eu.isas.searchgui.cmd.FastaCLI -in {} -decoy'.format(
                self.searchgui_path,
                self.fasta_db
                )
        self.db_path = os.path.join(self.db_cache, self.fasta_db.replace('.fasta', '_concatenated_target_decoy.fasta'))
        result = subprocess.run(cmd.split(), capture_output=True)
        print(result.stdout.decode())
        print(result.stderr.decode())

    def get_search_cmd(self, **kwargs):
        print(self.searchgui_path)
        cmd  = 'java -Xmx27G -cp {} eu.isas.searchgui.cmd.SearchCLI '.format(self.searchgui_path)
        cmd += '-spectrum_files {} '.format(self.mgf_path)
        cmd += '-output_folder {} '.format(self.out_dir)
        cmd += '-id_params {} '.format(self.id_params_path)
        for k,v in self.search_engines.items():
            cmd += '-{} {} '.format(k, v)
        cmd += '-protein_fdr {} '.format(self.protein_fdr)
        return cmd

class PeptideShaker:
    def __init__(self, searchgui, peptideshaker_version, sample_name = 'test', replicate   = 1):
        self.searchgui = searchgui
        self.tmp_dir = searchgui.tmp_dir
        self.out_dir = searchgui.out_dir
        self.peptideshaker_path = opj(self.searchgui.compomics_path, peptideshaker_version, f'{peptideshaker_version}.jar')
        self.sample_name = sample_name
        self.replicate   = replicate

        cmd = 'java -cp {} eu.isas.peptideshaker.cmd.PathSettingsCLI -temp_folder {} -ptm_configuration {}'.format(
                self.peptideshaker_path,
                self.tmp_dir,
                os.path.join(self.tmp_dir, 'ptmFactory-4.12.14.json')
                )
        result = subprocess.run(cmd, capture_output=True, shell=True)
        print(result.stdout.decode())
        print(result.stderr.decode())

    def run(self, ):
        cmd = 'java -Xmx27G -cp {} eu.isas.peptideshaker.cmd.PeptideShakerCLI '.format(self.peptideshaker_path)
        cmd += '-temp_folder {} '.format(self.tmp_dir)
        cmd += '-experiment {} '.format(self.searchgui.exp_name)
        cmd += '-sample {} '.format(self.sample_name)
        cmd += '-replicate {} '.format(self.replicate)
        cmd += '-id_params {} '.format(self.searchgui.id_params_path)
        cmd += '-identification_files {}/searchgui_out.zip '.format(self.searchgui.out_dir)
        cmd += '-spectrum_files {} '.format(self.searchgui.mgf_path)
        cmd += '-out {}/{}.cpsx '.format(self.out_dir, self.searchgui.exp_name)
        cmd += '-zip {}/{}.cpsx.zip '.format(self.out_dir, self.searchgui.exp_name)
        result = subprocess.run(cmd, capture_output=True, shell=True)
        print(result.stdout.decode())
        print(result.stderr.decode())

    def generate_reports(self, ):
        self.out_reports_dir = os.path.join(self.out_dir, 'peptideshaker_reports')
        if not os.path.exists(self.out_reports_dir):
            os.mkdir(self.out_reports_dir)
        reports_n = ', '.join([str(x) for x in range(9)])
        cmd  = 'java -cp {} eu.isas.peptideshaker.cmd.ReportCLI '.format(self.peptideshaker_path)
        cmd += '-in {}/{}.cpsx.zip '.format(self.out_dir, self.searchgui.exp_name)
        cmd += '-out_reports {} '.format(self.out_reports_dir)
        cmd += '-reports "{}"'.format(reports_n)
        result = subprocess.run(cmd, capture_output=True, shell=True)
        print(result.stdout.decode())
        print(result.stderr.decode())
