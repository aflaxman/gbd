from django.shortcuts import render_to_response, get_object_or_404
from django.contrib.auth.decorators import login_required
from django.http import *
from django.core.urlresolvers import reverse
from django.utils.translation import ugettext as _
from django import forms

from StringIO import StringIO

import gbd.fields
import gbd.view_utils as view_utils
from gbd.unicode_csv_reader import unicode_csv_reader
import dismod3

from models import *
from gbd.dismod3.utils import clean
from gbd.dismod3.settings import JOB_LOG_DIR, JOB_WORKING_DIR, SERVER_LOAD_STATUS_HOST, SERVER_LOAD_STATUS_PORT, SERVER_LOAD_STATUS_SIZE
from gbd.dismod3.table import population_by_region_year_sex
from gbd.dismod3.neg_binom_model import countries_for
import fcntl


class NewDataForm(forms.Form):
    file  = forms.FileField()
    required_data_fields = ['GBD Cause', 'Region', 'Parameter', 'Sex', 'Country ISO3 Code',
                            'Age Start', 'Age End', 'Year Start', 'Year End',
                            'Parameter Value', 'Units', ]

    tab_separated_values = \
        forms.CharField(required=False,
                        widget=forms.Textarea(attrs={'rows':20, 'cols':110, 'wrap': 'off'}),
                        help_text=_('See <a href="/public/file_formats.html">file format specification</a> for details, and consider using the <b>new</b> <a href="/public/datachecker.jnlp">Data Checker</a> app to clean your data.'))
    file  = forms.FileField(required=False)
    
    def clean_tab_separated_values(self):
        tab_separated_values = self.cleaned_data['tab_separated_values']
        if not tab_separated_values:
            if not self.files.has_key('file'):
                raise forms.ValidationError(_('TSV field and file field cannot both be blank'))
            else:
                return tab_separated_values
        lines = unicode_csv_reader(StringIO(tab_separated_values), dialect='excel-tab')
        return self.validate(lines)

    def clean_file(self):
        if self.file:
            file_data = self.file.read()
            lines = unicode_csv_reader(StringIO(file_data), dialect='excel-tab')
            return self.validate(lines)

    def validate(self, lines):
        """
Required data fields:
--------------------------------------------------------------------------------
Name                               Type    Limit
--------------------------------------------------------------------------------
GBD Cause                          str     one of the GBD causes
Region                             str     one of the GBD regions
Parameter                          str     standardize_data_type
Sex                                str     standardize_sex
Country ISO3 Code                  str     an ISO3 code in the region (or blank to apply to all countries in region)
Age Start                          int     [0, 100], <= Age End
Age End                            int     [0, 100], >= Age Start
Year Start                         int     [1950, 2010], <= Year End
Year End                           int     [1950, 2010], >= Year Start
Parameter Value                    float   >= 0
Units                              float   >= 1

Recommended data fields:
--------------------------------------------------------------------------------
Name                               Type    Limit
--------------------------------------------------------------------------------
Study ID                           empty or int     >= 0
Sequela                            empty or str     one of the GBD sequela codes
Case Definition                    empty or str     none
Coverage                           empty or float   [0,1]
Effective Sample Size*             empty or int     > 0, <= Total Study Size N
Lower CI*                          empty or float   >= 0 <= Parameter Value
Upper CI*                          empty or float   > Parameter Value
Standard Error*                    empty or float   > 0
Total Study Size N                 empty or int     > 0
Design Factor                      empty or float   >= 1
Citation                           empty or str     none
Urbanicity                         empty or float   [0, 1]
Ignore                             empty or int     [0, 1]

Optional data fields:
No checks

* Either of Effective Sample Size, Lower CI and Upper CI, or Standard Error must be given.
        """
        col_names = [clean(col) for col in lines.next()]

        # check that required fields appear
        for field in NewDataForm.required_data_fields:
            if not clean(field) in col_names:
                raise forms.ValidationError(_('Column "%s" is missing') % field)

        data_list = []
        for ii, cells in enumerate(lines):
            # skip blank lines
            if sum([cell == '' for cell in cells]) == len(cells):
                continue
            
            # ensure that something appears for each column
            if len(cells) != len(col_names):
                raise forms.ValidationError(
                    _('Error loading row %d:  found %d fields (expected %d))')
                    % (ii+2, len(cells), len(col_names)))

            # make an associative array from the row data
            data = {}
            for key, val in zip(col_names, cells):
                data[clean(key)] = val.strip()
            data['_row'] = ii+2

            data_list.append(data)

        # ensure that certain cells are the right format
        error_str = _('Row %d:  could not understand entry for %s')
        gbd_cause = ''

        for r in data_list:
            # check required data fields
            try:
                r['gbd_cause'] = str(r['gbd_cause'])
            except ValueError:
                raise forms.ValidationError(error_str % (r['_row'], 'GBD Cause'))
            if gbd_cause == '':
                gbd_cause = r['gbd_cause']
            else:
                if gbd_cause != r['gbd_cause']:
                    raise forms.ValidationError(error_str % (r['_row'], 'GBD Cause (all GBD Causes must be the same)'))

            try:
                r['region'] = str(r['region'])
            except ValueError:
                raise forms.ValidationError(error_str % (r['_row'], 'Region'))
            if not clean(r['region']) in [clean(region) for region in dismod3.gbd_regions] + ['all']:
                raise forms.ValidationError(error_str % (r['_row'], 'Region'))

            try:
                r['parameter'] = gbd.fields.standardize_data_type[r['parameter']]
            except KeyError:
                raise forms.ValidationError(error_str % (r['_row'], 'Parameter'))

            try:
                r['sex'] = gbd.fields.standardize_sex[r['sex']]
            except KeyError:
                raise forms.ValidationError(error_str % (r['_row'], 'Sex'))

            try:
                r['country_iso3_code'] = str(r['country_iso3_code'])
            except ValueError:
                raise forms.ValidationError(error_str % (r['_row'], 'Country ISO3 Code'))
            if r['region'] != 'all':
                if not r['country_iso3_code'] in countries_for[clean(r['region'])] + ['']:
                    raise forms.ValidationError(error_str % (r['_row'], 'Country ISO3 Code (%s is not in %s)' % (r['country_iso3_code'], r['region'])))
            elif r['country_iso3_code'] != 'all':
                raise forms.ValidationError(error_str % (r['_row'], 'Country ISO3 Code (%s must be "all" if region is "all")' % r['country_iso3_code']))

            try:
                r['age_start'] = int(r['age_start'])
            except ValueError:
                raise forms.ValidationError(error_str % (r['_row'], 'Age Start'))
            if r['age_start'] < 0 or r['age_start'] > 100:
                raise forms.ValidationError(error_str % (r['_row'], 'Age Start (must be in range [0, 100])'))

            try:
                r['age_end'] = int(r['age_end'])
            except ValueError:
                raise forms.ValidationError(error_str % (r['_row'], 'Age End'))
            if r['age_end'] < 0 or r['age_end'] > 100:
                raise forms.ValidationError(error_str % (r['_row'], 'Age End (must be in range [0, 100])'))

            if r['age_start'] > r['age_end']:
                raise forms.ValidationError(error_str % (r['_row'], 'Age Start (must be greater than Age End)'))

            try:
                r['year_start'] = int(r['year_start'])
            except ValueError:
                raise forms.ValidationError(error_str % (r['_row'], 'Year Start'))
            if r['year_start'] < 1950 or r['year_start'] > 2010:
                raise forms.ValidationError(error_str % (r['_row'], 'Year Start (must be in range [1950, 2010])'))

            try:
                r['year_end'] = int(r['year_end'])
            except ValueError:
                raise forms.ValidationError(error_str % (r['_row'], 'Year End'))
            if r['year_end'] < 1950 or r['year_end'] > 2010:
                raise forms.ValidationError(error_str % (r['_row'], 'Year End (must be in range [1950, 2010])'))
   
            if r['year_start'] > r['year_end']:
                raise forms.ValidationError(error_str % (r['_row'], 'Year Start (must be greater than Year End)'))

            units = 0
            try:
                units = float(r['units'].replace(',', '').replace('per ', ''))
            except ValueError:
                raise forms.ValidationError(error_str % (r['_row'], 'Units'))
            if units < 1:
                raise forms.ValidationError(error_str % (r['_row'], 'Units (must be greater than 1)'))

            try:
                r['parameter_value'] = float(r['parameter_value'])
            except ValueError:
                raise forms.ValidationError(error_str % (r['_row'], 'Parameter Value'))
            if r['parameter_value'] < 0:
                raise forms.ValidationError(error_str % (r['_row'], 'Parameter Value (must be greater than 0)'))
            value = r['parameter_value'] / units
            param = r['parameter']
            if param == 'prevalence data' and value > 1:
                raise forms.ValidationError(error_str % (r['_row'], 'Parameter Value of prevalence (must not be greater than 1)'))
            if param == 'duration data' and value > 100:
                raise forms.ValidationError(error_str % (r['_row'], 'Parameter Value of duration (must not be greater than 100)'))
            if param == 'relative-risk data' and value < 1:
                raise forms.ValidationError(error_str % (r['_row'], 'Parameter Value of relative-risk (must not be smaller than 1)'))
            if param == 'smr data' and value < 1:
                raise forms.ValidationError(error_str % (r['_row'], 'Parameter Value of smr (must not be smaller than 1)'))

            # check recommended data fields
            if 'study_id' in col_names and r['study_id'] != '':
                try:
                    r['study_id'] = int(r['study_id'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Study ID'))
                if r['study_id'] < 0:
                    raise forms.ValidationError(error_str % (r['_row'], 'Study ID (must be greater than 0)'))

            #if 'sequela' in col_names and r['sequela'] != '':
            #    try:
            #        r['sequela'] = str(r['sequela'])
            #    except ValueError:
            #        raise forms.ValidationError(error_str % (r['_row'], 'Sequela'))

            #if 'case_definition' in col_names and r['case_definition'] != '':
            #    try:
            #        r['case_definition'] = str(r['case_definition'])
            #    except ValueError:
            #        raise forms.ValidationError(error_str % (r['_row'], 'Case Definition'))

            if 'coverage' in col_names and r['coverage'] != '':
                try:
                    r['coverage'] = float(r['coverage'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Coverage'))
                if r['coverage'] < 0 or r['coverage'] > 1:
                    raise forms.ValidationError(error_str % (r['_row'], 'Coverage (must be in range [0, 1])'))

            effective_sample_size = 'effective_sample_size' in col_names and r['effective_sample_size'] != ''
            lower_ci = 'lower_ci' in col_names and r['lower_ci'] != ''
            upper_ci = 'upper_ci' in col_names and r['upper_ci'] != ''
            standard_error = 'standard_error' in col_names and r['standard_error'] != ''

            if not (effective_sample_size or (lower_ci and upper_ci) or standard_error):
                raise forms.ValidationError(error_str % (r['_row'], 'Either Effective Sample Size or both Lower CI and Upper CI or Standard Error must be given'))

            if effective_sample_size:
                try:
                    r['effective_sample_size'] = int(r['effective_sample_size'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Effective Sample Size'))
                if r['effective_sample_size'] <= 0:
                    raise forms.ValidationError(error_str % (r['_row'], 'Effective Sample Size (must be greater than 0)'))

            if lower_ci:
                try:
                    r['lower_ci'] = float(r['lower_ci'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Lower CI'))
                if r['lower_ci'] < 0 or r['lower_ci'] > r['parameter_value']:
                    raise forms.ValidationError(error_str % (r['_row'], 'Lower CI (must be less than parameter value)'))

            if upper_ci:
                try:
                    r['upper_ci'] = float(r['upper_ci'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Upper CI'))
                if r['upper_ci'] <= r['parameter_value']:
                    raise forms.ValidationError(error_str % (r['_row'], 'Upper CI (must be greater than Parameter Value)'))

            if standard_error:
                try:
                    r['standard_error'] = float(r['standard_error'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Standard Error'))
                if r['standard_error'] <= 0 and r['standard_error'] != -99:
                    raise forms.ValidationError(error_str % (r['_row'], 'Standard Error (must be greater than 0 or -99 for missing)'))

            if 'total_study_size_n' in col_names and r['total_study_size_n'] != '':
                try:
                    r['total_study_size_n'] = int(r['total_study_size_n'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Total Study Size N'))
                if r['total_study_size_n'] <= 0:
                    raise forms.ValidationError(error_str % (r['_row'], 'Total Study Size N (must be greater than 0)'))

            if 'total_study_size_n' in col_names and 'effective_sample_size' in col_names and r['effective_sample_size'] != '' and r['total_study_size_n'] != '':
                if r['effective_sample_size'] > r['total_study_size_n']:
                    raise forms.ValidationError(error_str % (r['_row'], 'Effective Sample Size (must be at most Total Study Size N)'))

            if 'design_factor' in col_names and r['design_factor'] != '':
                try:
                    r['design_factor'] = float(r['design_factor'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Design Factor'))
                if r['design_factor'] < 1:
                    raise forms.ValidationError(error_str % (r['_row'], 'Design Factor (must be greater than 1)'))

            #if 'citation' in col_names and r['citation'] != '':
            #    try:
            #        r['citation'] = str(r['citation'])
            #    except ValueError:
            #        raise forms.ValidationError(error_str % (r['_row'], 'Citation'))

            if 'urbanicity' in col_names and r['urbanicity'] != '':
                try:
                    r['urbanicity'] = float(r['urbanicity'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Urbanicity'))
                if r['urbanicity'] < 0 or r['urbanicity'] > 1:
                    raise forms.ValidationError(error_str % (r['_row'], 'Urbanicity (must be in range [0, 1])'))

            if 'ignore' in col_names and r['ignore'] != '':
                try:
                    r['ignore'] = int(r['ignore'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Ignore'))
                if r['ignore'] < 0 or r['ignore'] > 1:
                    raise forms.ValidationError(error_str % (r['_row'], 'Ignore (must be 0 or 1)'))

        return data_list

    
class NewDiseaseModelForm(forms.Form):
    model_json = \
        forms.CharField(required=True,
                        widget=forms.Textarea(attrs={'rows':20, 'cols':80, 'wrap': 'off'}),
                        help_text=_('See <a href="/public/dismod_data_json.html">dismod json specification</a> for details.'))
    def clean_model_json(self):
        model_json = self.cleaned_data['model_json']
        try:
            model_dict = json.loads(model_json)
        except ValueError:
            raise forms.ValidationError('JSON object could not be decoded')
        if not model_dict.get('params'):
            raise forms.ValidationError('missing params')
        if not model_dict.has_key('id'):
            raise forms.ValidationError('missing model id' % key)

        # store the model dict for future use
        self.cleaned_data['model_dict'] = model_dict
        return model_json
