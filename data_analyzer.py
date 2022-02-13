from scipy.stats import gamma
import datetime as dt
import geopandas as gpd
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
import sys

# Constant parameters
PLACES = ['Philippines', 'Baguio City', 'NCR', 'Philippines ex NCR']
INPUT_FILES = []
INPUT_PATH = '/mnt/c/Users/ralph/Downloads/'
DATE = dt.date.today()
date_str = DATE.strftime('%Y%m%d')
regex = r'^DOH.*{}.*Testing.*csv$|^DOH.*{}.*Case.*csv$'.format(date_str, date_str)
for file in os.listdir(INPUT_PATH):
    if re.search(regex, file):
        INPUT_FILES.append(os.path.join(INPUT_PATH, file))
CASE_INFO, TESTING_DATA = INPUT_FILES
config = {
    'user': 'rap',
    'password': 'paper',
    'host': os.environ.get("WSL_HOST_IP",""),
    'database': 'doh_data',
    'allow_local_infile': True
}

# Database tables
TABLES = {}
TABLES['testing_data'] = (
    "CREATE TABLE `testing_data` ("
    "   `name` VARCHAR(255) NOT NULL,"
    "   `date` CHAR(25) NOT NULL,"
    "   `turnaround` INT,"
    "   `samples` INT,"
    "   `tested` INT,"
    "   `positive` INT,"
    "   `negative` INT,"
    "   `equivocal` INT,"
    "   `invalid` INT,"
    "   `remain` INT,"
    "   `backlog` INT,"
    "   `cum_samples` INT,"
    "   `cum_tested` INT,"
    "   `cum_pos` INT,"
    "   `cum_neg` INT,"
    "   `pct_pos` INT,"
    "   `pct_neg` INT,"
    "   `validation` VARCHAR(255),"
    "   INDEX(`name`)"
    ")")

TABLES['case_info'] = (
    "CREATE TABLE `case_info` ("
    "   `case_code` VARCHAR(255) NOT NULL,"
    "   `age` INT,"
    "   `age_grp` CHAR(25),"
    "   `sex` CHAR(10),"
    "   `date_specimen` CHAR(25),"
    "   `date_release` CHAR(25),"
    "   `date_confirmed` CHAR(25) NOT NULL,"
    "   `date_died` CHAR(25),"
    "   `date_recovered` CHAR(25),"
    "   `removal_type` CHAR(25),"
    "   `admitted` CHAR(5),"
    "   `region` VARCHAR(255),"
    "   `PROVINCE` VARCHAR(255),"
    "   `city` VARCHAR(255),"
    "   `city_code` VARCHAR(255),"
    "   `brgy` VARCHAR(255),"
    "   `brgy_code` VARCHAR(255),"
    "   `status` CHAR(10),"
    "   `quarantine` CHAR(5),"
    "   `date_onset` CHAR(25),"
    "   `pregnant` CHAR(5),"
    "   `validation` VARCHAR(255),"
    "   INDEX(`region`, `city`, `province`)"
    ")")

# Drop table commands
commands = {}
commands['testing_data'] = 'DROP TABLE IF EXISTS testing_data'
commands['case_info'] = 'DROP TABLE IF EXISTS case_info'

# Imports data
imports = {}
imports['TESTING_DATA'] = (
    "LOAD DATA LOCAL INFILE '{}' "
    "INTO TABLE testing_data "
    "COLUMNS TERMINATED BY ',' "
    "OPTIONALLY ENCLOSED BY '{}' "
    "ESCAPED BY '{}' "
    "LINES TERMINATED BY '\n' IGNORE 1 ROWS"
    ).format(TESTING_DATA, '"', '"')

imports['CASE_INFO'] = (
    "LOAD DATA LOCAL INFILE '{}' "
    "INTO TABLE case_info "
    "COLUMNS TERMINATED BY ',' "
    "OPTIONALLY ENCLOSED BY '{}' "
    "ESCAPED BY '{}' "
    "LINES TERMINATED BY '\n' IGNORE 1 ROWS"
    ).format(CASE_INFO, '"', '"')

# Math constants
mu_R = 2.6
std_R = 2.0
mu_SI = 4.8
std_SI = 2.3
A = (mu_R / std_R) ** 2
B = (std_R**2) / mu_R
alpha = (mu_SI / std_SI)**2
beta = (std_SI**2) / mu_SI
r = np.arange(1,8)
w = gamma.cdf(r, a=alpha, scale=beta)

# NCR and Baguio City list of Laboratories
NCR_LAB = ('Asian Hospital and Medical Center', 'Chinese General Hospital',
           'Detoxicare Molecular Diagnostics Laboratory',
           'Dr. Jose N. Rodriguez Memorial Hospital and Sanitarium (TALA) GeneXpert Laboratory',
           'Lung Center of the Philippines (LCP)',
           'Lung Center of the Philippines GeneXpert Laboratory',
           'Makati Medical Center (MMC)', 'Marikina Molecular Diagnostics Laboratory (MMDL)',
           'UP Philippine Genome Center', 'Philippine Red Cross - Port Area',
           'Philippine Red Cross Logistics & Multipurpose Center', 'Philippine Red Cross (PRC)',
           'Research Institute for Tropical Medicine (RITM)', 'San Lazaro Hospital (SLH)',
           'Singapore Diagnostics',  "St. Luke's Medical Center - BGC (SLMC-BGC)",
           "St. Luke's Medical Center - BGC (HB) GeneXpert Laboratory",
           "St. Luke's Medical Center - Quezon City (SLMC-QC)",
           "St. Luke's Medical Center - Quezon City GeneXpert Laboratory", "The Medical City (TMC)",
           'UP National Institutes of Health (UP-NIH)', 'UP-PGH Molecular Laboratory',
           'National Kidney and Transplant Institute',
           'National Kidney and Transplant Institute GeneXpert Laboratory', 'PNP Crime Laboratory',
           'Safeguard DNA Diagnostics, Inc', 'Tondo Medical Center GeneXpert Laboratory',
           'Tropical Disease Foundation', 'Hi-Precision Diagnostics (QC)',
           'San Miguel Foundation Testing Laboratory',
           'Sta. Ana Hospital - Closed System Molecular Laboratory (GeneXpert)',
           'Sta. Ana Hospital - Closed System Molecular Laboratory (RT PCR)',
           'Fe del Mundo Medical center', 'Health Delivery Systems',
           'Dr. Jose N. Rodriguez Memorial Hospital and Sanitarium (TALA) RT PCR',
           'Veteran Memorial Medical Center', 'Philippine Heart Center GeneXpert Laboratory',
           'Philippine Airport Diagnostic Laboratory', 'De Los Santos Medical Center',
           'Victoriano Luna - AFRIMS', "Philippine Children's Medical Center",
           "University of Perpetual Help DALTA Medical Center, Inc.",
           'AL Molecular Diagnostic Laboratory', "Amosup Seamen's Hospital",
           "Valenzuela Hope Molecular Laboratory",
           'Las Pinas General Hospital and Satellite Trauma Center',
           'Las Pinas General Hospital and Satellite Trauma Center GeneXpert Laboratory',
           'Kaiser Medical Center Inc.',
           'Taguig City Molecular Laboratory', 'First Aide Diagnostic Center',
           'New World Diagnostic Premium Medical Branch', 'BioPath Clinical Diagnostic, Inc.',
           "The Lord's Grace Medical and Industrial Clinic", "Manila Doctors Hospital",
           'The Premier Molecular Diagnostics', 'Cardinal Santos Medical Center',
           'Army General Hospital Molecular Laboratory',
           'University of the East Ramon Magsaysay Memorial Medical Center',
           'Health Metrics', 'IOM - Manila Health Centre Laboratory - GeneXpert',
           'PNP General Hospital', 'QR Medical Laboratories Inc. (VitaCare)',
           'South Super Hi Way Molecular Diagnostic Lab',
           'Manila Healthtek Inc.', 'BioPath Clinical Diagnostics, Inc. (E. Rodriguez)',
           'Rizal Medical Center GeneXpert Laboratory', 'Kairos Diagnostics Laboratory',
           'Amang Rodriguez Memorial Center GeneXpert Laboratory',
           'Hero Laboratories', 'Lifecore Biointegrative Inc.',
           'BioPath Clinical Diagnostics, Inc. GeneXpert Laboratory',
           'Pasig City Molecular Laboratory', 'Amang Rodriguez Memoral Medical Center (RT PCR)',
           'BioPath Clinical Diagnostics, Inc. (E. Rodriguez) GeneXpert Laboratory',
           'Manila Diagnostic Center for OFW Inc.', 'Supercare Medical Services, Inc.',
           'Philippine General Hospital GeneXpert Laboratory', 'Quezon City Molecular Diagnostics Laboratory',
           'JT Cenica Medical Health System', 'Caloocan City North Medical Center', 'A Star Laboratories',
           'Jose R. Reyes Memorial Medical Center GeneXpert Laboratory', 'Ospital ng Muntinlupa',
           'Chinese General Hospital - GeneXpert Laboratory',
           'Philippine Medical Diagnostic & Laboratory Center Corporation',
           'The Meidcal City (TMC) - GeneXpert Laboratory', 'Marilao Medical and Diagnosc Clinic, Inc. - Pasay City',
           'New World Diagnostic Premium Medical Branch - GeneXpert Laboratory',
           'El Roi Molecular Diagnostic Laboratory', 'Quirino Memorial Medical Center',
           'Ateneo Molecular Pathology Laboratory', 'Victoriano Luna - AFRIMS GeneXpert Laboratory',
           'Truehealth Diagnostics Center Corp.', 'PSG Station Hospital Molecular Laboratory',
           'Ospital ng Parañaque District II', 'Olayn Medical Laboratory, Inc. - GeneXpert Laboratory',
           'Olayn Medical Laboratory, Inc.', 'Multi-Lifecare Foundation Molecular Laboratory',
           'Metro North Medical Center and Hospital Inc.',
           'Medical Center Parañaque, Inc.', 'Medical Center of Taguig City, Inc.',
           'Intramuros Molecular Laboratory', 'Interpharma Solutions Philippines, Inc.',
           'Hi-Precision Diagnostics Cartridge-Based PCR Laboratory (QC)', 'Heartland Medical and Diagnostic Center',
           'Healthway Medical', 'F.Y Manalo Medical Foundation, Inc. (formerly New Era Hospital)',
           'Ermita Molecular Diagnostic Laboratory, Inc.', 'East Avenue Medical Center',
           'Cardinal Santos Medical Center Cartridge-based PCR', 'Best Diagnostic Corporation',
           'Allied Care Expert - Valenzuela', 'Kairos Diagnostics Laboratory Cartridge-based PCR',
           'AL Molecular Diagnostic Laboratory - GeneXpert', 'Marikina Valley Medical Center')
BAGUIO_LAB = ('Baguio General Hospital and Medical Center', 'Parkway Medical and Diagnostic Center')

# Directory where to store program outputs
DIR = '/mnt/g/My Drive/Desktop_Dump/'
FILE_DIR = os.path.join(DIR, 'COVID_Analysis_PosRate_RepNum')
if not os.path.exists(FILE_DIR):
    os.mkdir(FILE_DIR)

def reproduction_number(I):
    I_s = [np.sum(I[0:7-i]) for i in range(0,7)]
    I_s = np.array(I_s)
    Delta = np.dot(w,I_s)
    num = (A + np.sum(I[1:]))
    R = num / (1/B + Delta)
    CV = 1 / math.sqrt(num)
    CI = 1.96 * CV * R
    return R, CI

def analyzer(location, sql_con, sql_cursor):
    # Date parameters for R
    start_date = dt.date(2020,3,16)
    start_onset = dt.date(2020,3,8)
    end_date = DATE
    tau = dt.timedelta(days=7)

    # Date range for plots
    datemin = dt.date(2020,3,1)
    datemax = dt.date(2021,12,1)
    date_range = [datemin, datemax]
    fill_dates = pd.date_range(start=datemin, end=datemax)

    # File directory of output files
    OUTPUT = os.path.join(FILE_DIR, 'doh_data_{}.csv'.format(location))
    IMAGE = os.path.join(FILE_DIR, 'pos_rate_{}.png'.format(location))
    EPI_CURVE = os.path.join(FILE_DIR, 'epi_curve_{}.png'.format(location))
    PROVINCIAL_CASES = os.path.join(FILE_DIR, 'provincial_case_count.csv')
    NCR_CASES = os.path.join(FILE_DIR, 'ncr_case_count.csv')
    CUMMULATIVE_CASES = os.path.join(FILE_DIR, 'cummulative_cases.csv')
    BENGUET_DATA = os.path.join(FILE_DIR, 'benguet_cummulative.csv')

    # Set laboratory list depending on location
    if location == 'Philippines':
        LAB = []
        file_url = os.path.join(DIR, 'lab_list.txt')
        with open(file_url,'r+') as lab_file:
            lab_list = lab_file.read().split('\n')
            lab_list.remove('')
            #for lab in sorted(lab_list):
            #    print(lab)
            print('Length old lab list: {}'.format(len(lab_list)))
            query = "SELECT DISTINCT name FROM testing_data ORDER BY name"
            sql_cursor.execute(query)
            count = 0
            for lab in sql_cursor:
                LAB.append(lab[0])
                if lab[0] not in lab_list:
                    lab_file.write('{}\n'.format(lab[0]))
                    print('New laboratory: {}'.format(lab[0]))
                    count += 1
            print('Number of new laboratories: {}'.format(count))
        case_query = (
            "SELECT date_onset AS date, date_confirmed AS dconf, date_specimen AS dspec, "
            "date_release AS drel FROM case_info"
            )
        death_query = "SELECT date_confirmed AS date, date_died AS ddied FROM case_info WHERE removal_type='DIED'"
        #for lab in sorted(LAB):
        #    print(lab)
        #sys.exit(1)

        # Extract cummulative cases, deaths and recoveries
        total_query = (
            "SELECT DISTINCT PROVINCE, COUNT(date_confirmed) AS total FROM case_info GROUP BY PROVINCE"
        )
        recoveries_query = (
            "SELECT DISTINCT PROVINCE, COUNT(removal_type) AS recoveries FROM case_info "
            "WHERE removal_type='RECOVERED' GROUP BY PROVINCE"
        )
        deaths_query = (
            "SELECT DISTINCT PROVINCE, COUNT(removal_type) AS deaths FROM case_info "
            "WHERE removal_type='DIED' GROUP BY PROVINCE"
        )

        total = pd.read_sql(total_query, sql_con, index_col=['PROVINCE'])
        recoveries = pd.read_sql(recoveries_query, sql_con, index_col=['PROVINCE'])
        deaths = pd.read_sql(deaths_query, sql_con, index_col=['PROVINCE'])
        cummulative = pd.concat([total, recoveries, deaths], axis=1)
        cummulative['deaths'].fillna(0, inplace=True)
        cummulative.loc['MAGUINDANAO','total'] += cummulative.loc['COTABATO CITY (NOT A PROVINCE)', 'total']
        cummulative.loc['MAGUINDANAO','recoveries'] += cummulative.loc['COTABATO CITY (NOT A PROVINCE)', 'recoveries']
        cummulative.loc['MAGUINDANAO','deaths'] += cummulative.loc['COTABATO CITY (NOT A PROVINCE)', 'deaths']
        cummulative.loc['BASILAN','total'] += cummulative.loc['CITY OF ISABELA (NOT A PROVINCE)', 'total']
        cummulative.loc['BASILAN','recoveries'] += cummulative.loc['CITY OF ISABELA (NOT A PROVINCE)', 'recoveries']
        cummulative.loc['BASILAN','deaths'] += cummulative.loc['CITY OF ISABELA (NOT A PROVINCE)', 'deaths']
        cummulative.loc['DAVAO DEL SUR','total'] += cummulative.loc['DAVAO OCCIDENTAL', 'total']
        cummulative.loc['DAVAO DEL SUR','recoveries'] += cummulative.loc['DAVAO OCCIDENTAL', 'recoveries']
        cummulative.loc['DAVAO DEL SUR','deaths'] += cummulative.loc['DAVAO OCCIDENTAL', 'deaths']
        cummulative['active'] = cummulative['total'] - cummulative['recoveries'] - cummulative['deaths']
        cummulative = cummulative.drop(['COTABATO CITY (NOT A PROVINCE)', 'CITY OF ISABELA (NOT A PROVINCE)',
                                      'DAVAO OCCIDENTAL'])
        cummulative.index = cummulative.index.str.lower()
        change = cummulative.index.tolist()
        idx = change.index('cotabato (north cotabato)')
        change[idx] = 'north cotabato'
        idx = change.index('ncr')
        change[idx] = 'metropolitan manila'
        idx = change.index('samar (western samar)')
        change[idx] = 'samar'
        idx = change.index('davao de oro')
        change[idx] = 'compostela valley'
        cummulative.index = change
        cummulative.loc['shariff kabunsuan'] = cummulative.loc['maguindanao',:]
        cummulative = cummulative.sort_index()
        population = pd.read_csv(os.path.join(FILE_DIR, 'PH_pop.csv'), index_col=['PROVINCE'])
        population.index = population.index.str.lower()
        cummulative = pd.concat([cummulative, population], axis=1)
        cummulative['act_per_pop'] = cummulative['active'] / cummulative['POPULATION'] * 100000.0
        cummulative = cummulative.drop([""])
        if os.path.exists(CUMMULATIVE_CASES):
            os.remove(CUMMULATIVE_CASES)
        cummulative.to_csv(CUMMULATIVE_CASES)

        # Extract current case count for each province
        prov_query = (
            "SELECT DISTINCT PROVINCE, COUNT(date_confirmed) AS COUNT FROM case_info "
            "WHERE date_confirmed = {} GROUP BY PROVINCE"
        ).format("'{}'".format(DATE.strftime('%Y-%m-%d')))
        prov_cases = pd.read_sql(prov_query, sql_con, index_col=['PROVINCE'])
        prov_cases.index = prov_cases.index.str.lower()
        prov_cases = prov_cases.drop([""])

        # Extract current case count for NCR
        NCR_query = (
            "SELECT DISTINCT city, COUNT(date_confirmed) AS COUNT FROM case_info "
            "WHERE region = 'NCR' AND date_confirmed = {} GROUP BY city"
        ).format("'{}'".format(DATE.strftime('%Y-%m-%d')))
        NCR_cases = pd.read_sql(NCR_query, sql_con, index_col=['city'])
        NCR_cases.index = NCR_cases.index.str.lower()
        #NCR_cases = NCR_cases.drop([""])

        # Write provincial case counts to csv
        if os.path.exists(PROVINCIAL_CASES):
            os.remove(PROVINCIAL_CASES)
        prov_cases.to_csv(PROVINCIAL_CASES)

        # Write NCR case counts to csv
        if os.path.exists(NCR_CASES):
            os.remove(NCR_CASES)
        NCR_cases.to_csv(NCR_CASES)
    elif location == 'Baguio City':
        LAB = sorted(BAGUIO_LAB)
        case_query = (
            "SELECT date_onset AS date, date_confirmed AS dconf, date_specimen AS dspec, "
            "date_release AS drel FROM case_info WHERE city = {}"
            ).format("'{}'".format(location))
        death_query = '''
            SELECT date_confirmed AS date, date_died AS ddied FROM case_info
            WHERE removal_type = 'DIED' AND city = {}'''.format("'{}'".format(location))
        start_date = dt.date(2020,4,12)
        start_onset = start_date

        # Extract cummulative data for Benguet
        total_query = (
            "SELECT city, COUNT(date_confirmed) AS total FROM case_info WHERE PROVINCE='BENGUET' GROUP BY city"
        )
        recoveries_query = (
            "SELECT city, COUNT(removal_type) AS recoveries FROM case_info "
            "WHERE PROVINCE='BENGUET' AND removal_type='RECOVERED' GROUP BY city"
        )
        deaths_query = (
            "SELECT city, COUNT(removal_type) AS deaths FROM case_info "
            "WHERE PROVINCE='BENGUET' AND removal_type='DIED' GROUP BY city"
        )

        total = pd.read_sql(total_query, sql_con, index_col=['city'])
        recoveries = pd.read_sql(recoveries_query, sql_con, index_col=['city'])
        deaths = pd.read_sql(deaths_query, sql_con, index_col=['city'])
        cummulative = pd.concat([total, recoveries, deaths], axis=1)
        cummulative['deaths'].fillna(0, inplace=True)
        cummulative['active'] = cummulative['total'] - cummulative['recoveries'] - cummulative['deaths']
        cummulative.index = cummulative.index.str.lower()
        cummulative = cummulative.drop([""])
        population = pd.read_csv(os.path.join(FILE_DIR, 'benguet_pop.csv'), index_col=['city'])
        population.index = population.index.str.lower()
        benguet_data = pd.concat([cummulative, population], axis=1)
        benguet_data['act_per_pop'] = benguet_data['active'] / benguet_data['pop'] * 100000.0
        benguet_data['NAME'] = benguet_data.index.str.capitalize()
        benguet_data.loc['baguio city', 'NAME'] = 'Baguio City'
        benguet_data.loc['la trinidad (capital)', 'NAME'] = 'La Trinidad'
        if os.path.exists(BENGUET_DATA):
            os.remove(BENGUET_DATA)
        benguet_data.to_csv(BENGUET_DATA)
    elif location == 'NCR':
        LAB = sorted(NCR_LAB)
        case_query = (
            "SELECT date_onset AS date, date_confirmed AS dconf, date_specimen AS dspec, "
            "date_release AS drel FROM case_info WHERE region = {}"
            ).format("'{}'".format(location))
        death_query = '''
            SELECT date_confirmed AS date, date_died AS ddied FROM case_info
            WHERE removal_type = 'DIED' AND region = {}'''.format("'{}'".format(location))
    elif location == 'Philippines ex NCR':
        LAB = []
        query = (
            "SELECT DISTINCT name FROM testing_data WHERE name NOT IN {} ORDER BY name"
            ).format(NCR_LAB)
        sql_cursor.execute(query)
        for lab in sql_cursor:
            LAB.append(lab[0])
        case_query = (
            "SELECT date_onset AS date, date_confirmed AS dconf, date_specimen AS dspec, "
            "date_release AS drel FROM case_info WHERE region != {}"
            ).format("'NCR'")
        death_query = '''
            SELECT date_confirmed AS date, date_died AS ddied FROM case_info
            WHERE removal_type = 'DIED' AND region != {}'''.format("'NCR'")
    print('Total laboratories {}: {}'.format(location, len(LAB)))

    # Extract pertinent data from database
    # Testing and death statistics
    query = (
        "SELECT date, SUM(tested) AS Tested, "
        "SUM(positive) AS Positive, "
        "SUM(positive) / SUM(tested) * 100.0 AS `Positivity Rate` "
        "FROM testing_data WHERE name IN {} "
        "GROUP BY date"
        ).format(tuple(LAB))
    deaths = pd.read_sql(death_query, sql_con, parse_dates=['date','ddied'])
    for i in deaths.index:
        if pd.isnull(deaths.at[i, 'ddied'].to_pydatetime()):
            deaths.at[i, 'ddied'] = deaths.at[i, 'date']
    deaths = deaths.groupby('ddied').size().to_frame(name='Deaths')
    deaths.index.names = ['date']
    testing = pd.read_sql(query, sql_con, parse_dates=['date'], index_col=['date'])
    testing['Positivity Rate'].fillna(0, inplace=True)
    # Daily case count
    raw = pd.read_sql(case_query, sql_con, parse_dates=['date','dconf','dspec','drel'])
    for i in raw.index:
        if pd.isnull(raw.at[i, 'date'].to_pydatetime()):
            if not pd.isnull(raw.at[i,'dspec'].to_pydatetime()):
                raw.at[i, 'date'] = raw.at[i,'dspec']
    cases = raw.groupby('dconf').size().to_frame(name='New Cases')
    cases.index.names = ['date']
    cases_onset = raw.groupby('date').size().to_frame(name='Cases Date Onset')

    # Process with pandas for visualization and analysis
    final = pd.concat([testing, cases, cases_onset, deaths], axis=1)
    final['New Cases'].fillna(0, inplace=True)
    final['Cases Date Onset'].fillna(0, inplace=True)
    final['Positivity Rate'].dropna(inplace=True)
    final['R'] = np.nan
    final['lwr95'] = np.nan
    final['upr95'] = np.nan
    final['R_onset'] = np.nan
    final['lwr95_onset'] = np.nan
    final['upr95_onset'] = np.nan
    final['Deaths'].fillna(0, inplace=True)
    sma = final['Positivity Rate'].rolling(window=7).mean()
    sma_num_tests = final['Tested'].rolling(window=7).mean()
    sma_epi_curve = final['New Cases'].rolling(window=7).mean()
    sma_deaths = final['Deaths'].rolling(window=7).mean()
    upper = 5.0*np.ones(fill_dates.size)
    lower = np.zeros(fill_dates.size)

    # Compute R(t) based on date of confirmation
    #target_range = pd.date_range(start=start_date, end=end_date)
    #for day in target_range:
    #    case_range = pd.date_range(start=day-tau, end=day)
    #    I = final.loc[case_range,'New Cases'].values
    #    rep_num, CI = reproduction_number(I)
    #    final.at[day, ['R','lwr95', 'upr95']] = rep_num, rep_num - CI, rep_num + CI

    # Compute R(t) based on date of symptom onset
    #target_range = pd.date_range(start=start_onset, end=end_date)
    #for day in target_range:
    #    case_range = pd.date_range(start=day-tau, end=day)
    #    I = final.loc[case_range,'Cases Date Onset'].values
    #    rep_num, CI = reproduction_number(I)
    #    final.at[day, ['R_onset','lwr95_onset', 'upr95_onset']] = rep_num, rep_num - CI, rep_num + CI

    # Plot positivity rate and daily number of tests
    #fig, (ax, ax1) = plt.subplots(2,1)
    fig, (ax, ax1, ax2) = plt.subplots(3,1)
    fig.subplots_adjust(hspace=0.5)
    ax.bar(final.index.values, final['New Cases'],
            label='Date of reporting')
    #ax.bar(final.index.values, final['Cases Date Onset'],
    #        label='Onset of symptoms')
    ax.plot(final.index.values, sma_epi_curve.values, label='7-day MA',
             color='red')
    ax.set_xlim(date_range)
    ax.legend(loc='upper left')
    ax.set_xlabel('Date')
    ax.set_ylabel('Number of reported cases')
    ax.set_title('COVID-19 Epidemic Curve ({})'.format(location))

    ax1.plot(final.index.values, final['Positivity Rate'].values,
            color='orange',label='Positivity Rate (%)')
    ax1.plot(final.index.values, sma.values, color='blue',label='7-day MA')
    ax1.fill_between(fill_dates,lower,upper,alpha=0.3,color='green', label='WHO safe zone')
    ax1.set_ylim(0,50)
    ax1.set_xlim(date_range)
    ax1.set_xlabel('Date')
    ax1.legend(loc='upper left')
    ax1.set_title('Daily Positivity Rate ({})'.format(location))

    ax2.bar(final.index.values, final['Tested'].values, label='Tested')
    ax2.bar(final.index.values, final['Positive'].values, label='Positive')
    ax2.plot(final.index.values, sma_num_tests.values, label='7-day MA', color='red')
    ax2.set_ylim(0,math.ceil(final['Tested'].max()))
    ax2.set_xlim(date_range)
    ax2.set_xlabel('Date')
    ax2.set_ylabel('Number of individuals')
    ax2.legend(loc='upper left')
    ax2.set_title('Daily Testing Numbers ({})'.format(location))

    # Save testing plot to image
    if os.path.exists(IMAGE):
        os.remove(IMAGE)
    plt.gcf().set_size_inches(16.5,11.75)
    plt.savefig(IMAGE,dpi=600)

    # Plot epi_curve data
    #fig1, (ax2, ax3, ax4) = plt.subplots(3,1)
    #fig1, (ax2, ax3) = plt.subplots(2,1)
    #fig1.subplots_adjust(hspace=0.5)
    #ax2.bar(final.index.values, final['New Cases'],
    #        label='Date of reporting', color='orange')
    #ax2.bar(final.index.values, final['Cases Date Onset'],
    #        label='Onset of symptoms')
    #ax2.plot(final.index.values, sma_epi_curve.values, label='7-day MA',
    #         color='red')
    #ax2.set_xlim(date_range)
    #ax2.legend(loc='upper left')
    #ax2.set_xlabel('Date')
    #ax2.set_ylabel('Number of reported cases')
    #ax2.set_title('COVID-19 Epidemic Curve ({})'.format(location))

    # Plot Death curve
    #ax3.bar(final.index.values, final['Deaths'].values, color='grey', label='Deaths')
    #ax3.plot(final.index.values, sma_deaths.values, color='black', label='7-day MA')
    #ax3.set_xlim(date_range)
    #ax3.legend(loc='upper left')
    #ax3.set_xlabel('Date')
    #ax3.set_ylabel('Number of reported deaths')
    #ax3.set_title('Death Curve ({})'.format(location))

    # Plot Rt curve (date of reporting)
    #ax3.plot(final.index.values, final['R'].values, label='R-value', color='red')
    #ax3.fill_between(final.index.values,final['lwr95'], final['upr95'],
    #                 alpha=0.4, label='95% CI')
    #ax3.plot(fill_dates.values, np.ones(fill_dates.values.shape), 'k--')
    #ax3.set_xlim(date_range)
    #minimum = final['lwr95'].min()
    #if final['lwr95'].min() < 0:
    #    minimum = 0.0
    #ax3.set_ylim(minimum, 3.0) #final['upr95'].max())
    #ax3.set_xlabel('Date of reporting')
    #ax3.set_ylabel('Reproduction Number, R')
    #ax3.legend(loc='upper right')
    #ax3.set_title('Reproduction Number Based on Date of Reporting ({})'.format(location))

    # Plot Rt curve (date of symptom onset)
    #ax4.plot(final.index.values, final['R_onset'].values, label='R-value', color='green')
    #ax4.fill_between(final.index.values,final['lwr95_onset'], final['upr95_onset'],
    #                 alpha=0.4, label='95% CI')
    #ax4.plot(fill_dates.values, np.ones(fill_dates.values.shape), 'k--')
    #ax4.set_xlim(date_range)
    #ax4.set_ylim(0.0, final['upr95_onset'].max())
    #ax4.set_xlabel('Date of symptom onset')
    #ax4.set_ylabel('Reproduction Number, R')
    #ax4.legend(loc='upper right')
    #ax4.set_title('Reproduction Number Based on Date of Symptom Onset ({})'.format(location))

    # Save epi_curve plot to image
    #if os.path.exists(EPI_CURVE):
    #    os.remove(EPI_CURVE)
    #plt.gcf().set_size_inches(16.5,11.75)
    #plt.gcf().set_size_inches(11.75,8.25)
    #plt.savefig(EPI_CURVE,dpi=600)

    # Write processed data to csv
    if os.path.exists(OUTPUT):
        os.remove(OUTPUT)
    final.to_csv(OUTPUT)

    #if location == 'Philippines':
    #    sys.exit(1)

def run_sql_command(commands, method, cursor):
    if method == 'DROP':
        msg = 'Deleting table {}: '
    elif method == 'CREATE':
        msg = 'Creating table {}: '
    else:
        msg = 'Importing {}: '
    for command in commands:
        print(msg.format(command), end='')
        cursor.execute(commands[command])
        print('OK')

def data_mapper():
    # INPUT FILES
    CUM_DATA = os.path.join(FILE_DIR, 'cummulative_cases.csv')
    BENGUET_CUM = os.path.join(FILE_DIR, 'benguet_cummulative.csv')

    # OUTPUT FILES
    ACT_POP = os.path.join(FILE_DIR, 'provincial_active_per_pop_{}.png'.format(DATE))
    BEN_ACT_POP = os.path.join(FILE_DIR, 'benguet_heatmap_{}.png'.format(DATE))

    # Load data
    map_data = gpd.read_file('/mnt/g/My Drive/GIS_DATA_PH/Provinces.shp', crs='EPSG:4326')
    benguet_map = gpd.read_file('/mnt/g/My Drive/GIS_DATA_PH/Benguet.shp', crs='EPSG:4326')
    cummulative = pd.read_csv(CUM_DATA)
    cummulative = cummulative.set_index(cummulative.columns[0])
    benguet_cum = pd.read_csv(BENGUET_CUM, index_col=['city'])

    # Data clean up
    map_data['PROVINCE'] = map_data['PROVINCE'].str.lower()
    map_data = map_data.dissolve(by='PROVINCE')
    benguet_map = benguet_map.dissolve(by='TOWN')
    benguet_map.index = benguet_map.index.str.lower()

    # Merge data, map and population
    cum_data = pd.concat([map_data, cummulative], axis=1)
    benguet_map = pd.concat([benguet_map, benguet_cum],axis=1)
    benguet_map['coords'] = benguet_map['geometry'].apply(lambda x: x.representative_point().coords[:])
    benguet_map['coords'] = [coords[0] for coords in benguet_map['coords']]

    # Plot map
    cum_data.plot(column='act_per_pop', cmap='hot_r', legend=True, edgecolor='gray', linewidth=0.4,
                  legend_kwds={'label': 'Active cases per 100,000 population'}, vmax=500, vmin=0)
    ax = plt.gca()
    ax.set_facecolor((0.8,0.8,0.8))
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    colorbar = ax.get_figure().get_axes()[1]
    colorbar.set_yticklabels([0.0, 100, 200, 300, 400, r'$\geq500$'])
    ax.set_title('{}'.format(DATE))
    plt.gcf().set_size_inches(6,8.25)
    plt.savefig(ACT_POP, dpi=600)

    benguet_map.plot(column='act_per_pop', cmap='hot_r', legend=True, edgecolor='black', linewidth=0.4,
                     legend_kwds={'label': 'Active cases per 100,000 population'}, vmax=500, vmin=0)
    ax = plt.gca()
    ax.set_facecolor((0.8,0.8,0.8))
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    colorbar = ax.get_figure().get_axes()[1]
    colorbar.set_yticklabels([0.0, 100, 200, 300, 400, r'$\geq500$'])
    for idx, row in benguet_map.iterrows():
        plt.annotate(s=row['NAME'], xy=row['coords'],horizontalalignment='center',fontsize=6)
    ax.set_title('{}'.format(DATE))
    plt.gcf().set_size_inches(6,8.25)
    plt.savefig(BEN_ACT_POP, dpi=600)
