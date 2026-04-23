from logging import getLogger
import requests
import webbrowser
import time
import os
import json
import sqlite3
from collections import defaultdict
import numpy as np

from calour.database import Database
from calour.uniref_experiment import UniRefExperiment

logger = getLogger(__name__)


class Uniref(Database):
    '''uniref calour interface
    '''
    def __init__(self, exp=None):
        '''Called every time a database interface is created (e.g. when creating a plot, etc.)
        Can put here obtaining the database api address, handshake, etc.

        Note: the methods variable tells calour what things can be done with this database. can include:
            'get' if database interface supports get_seq_annotation_strings()
                This means the database returns a list of strings describing the selcted feature.
                used for the interavtive heatmap display when clicking a feature.
            'annotate' if database interface supports add_annotation()
                This means the user can click the "annotate" button to add selected sequences to the database. The GUI is fully
                provided by the database module.
            'enrichment' if database interface supports get_feature_terms()
                This means the database supports statistical tests to identify interesting properties about a group of features.
                Used for example in the "erichment" button in the heatmap.
        In minimum, the database should support the 'get' method.

        Parameters
        ----------
        exp: calour.Experiment
            the calour Experiment this database is associated with.
            useful if you want to store data obtained from the database in the experiment to prevent multiple database queries
            about the same data. In that case, can store the data in exp.exp_metadata['__redbiom_XXX']
        '''
        super().__init__(exp=exp, database_name='uniref', methods=['get'])
        self.uniref_url = 'https://rest.uniprot.org/'
        self.quickgo_url = 'https://www.ebi.ac.uk/QuickGO/services'
        self.cache_db_path = os.path.expanduser(
            os.environ.get('UNIREF_CALOUR_CACHE_DB', '~/.cache/uniref_calour/uniref_cache.sqlite')
        )
        self._init_cache_db()
        logger.debug('uniref calour database interface initialized')

    def set_log_level(self, level):
        '''Set the log level for this database interface.

        Parameters
        ----------
        level : int or str
            The log level to set (e.g., logging.DEBUG, logging.INFO, 'DEBUG', 'INFO').
        '''
        logger.setLevel(level)

    def _init_cache_db(self):
        '''Initialize the local sqlite cache database (persistent across sessions).'''
        cache_dir = os.path.dirname(self.cache_db_path)
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)
        with sqlite3.connect(self.cache_db_path) as conn:
            conn.execute(
                '''
                CREATE TABLE IF NOT EXISTS uniref_cache (
                    uniref_id TEXT PRIMARY KEY,
                    go_terms_json TEXT NOT NULL,
                    names_json TEXT NOT NULL,
                    accessions_json TEXT NOT NULL,
                    organisms_json TEXT NOT NULL DEFAULT '[]',
                    length INTEGER NOT NULL,
                    updated_at INTEGER NOT NULL
                )
                '''
            )
            conn.execute(
                '''
                CREATE TABLE IF NOT EXISTS go_term_cache (
                    go_id TEXT PRIMARY KEY,
                    go_name TEXT,
                    updated_at INTEGER NOT NULL
                )
                '''
            )
            conn.commit()

    def _get_cached_go_term_names(self, go_ids: list[str]) -> dict[str, str | None]:
        '''Return cached GO term names for the given ids.'''
        if not go_ids:
            return {}

        placeholders = ','.join(['?'] * len(go_ids))
        with sqlite3.connect(self.cache_db_path) as conn:
            rows = conn.execute(
                f'''
                SELECT go_id, go_name
                FROM go_term_cache
                WHERE go_id IN ({placeholders})
                ''',
                tuple(go_ids),
            ).fetchall()

        return {row[0]: row[1] for row in rows}

    def _set_cached_go_term_names(self, go_term_names: dict[str, str | None]):
        '''Upsert GO term names into the persistent cache.'''
        if not go_term_names:
            return

        now = int(time.time())
        values = [(go_id, go_name, now) for go_id, go_name in go_term_names.items()]

        with sqlite3.connect(self.cache_db_path) as conn:
            conn.executemany(
                '''
                INSERT INTO go_term_cache (go_id, go_name, updated_at)
                VALUES (?, ?, ?)
                ON CONFLICT(go_id) DO UPDATE SET
                    go_name = excluded.go_name,
                    updated_at = excluded.updated_at
                ''',
                values,
            )
            conn.commit()

    def _get_cached_uniref_info(self, uniref_id: str):
        '''Return cached UniRef info for the id, or None if not present.'''
        with sqlite3.connect(self.cache_db_path) as conn:
            row = conn.execute(
                '''
                SELECT go_terms_json, names_json, accessions_json, organisms_json, length
                FROM uniref_cache
                WHERE uniref_id = ?
                ''',
                (uniref_id,),
            ).fetchone()

        if row is None:
            return None

        return {
            'go_terms': json.loads(row[0]),
            'names': json.loads(row[1]),
            'accessions': json.loads(row[2]),
            'organisms': json.loads(row[3]),
            'length': int(row[4]),
        }

    def _set_cached_uniref_info(self, uniref_id: str, data: dict):
        '''Upsert UniRef info into the local cache.'''
        go_terms = data.get('go_terms', [])
        names = data.get('names', [])
        accessions = data.get('accessions', [])
        organisms = data.get('organisms', [])
        seq_len = int(data.get('length', 0) or 0)

        with sqlite3.connect(self.cache_db_path) as conn:
            conn.execute(
                '''
                INSERT INTO uniref_cache (
                    uniref_id,
                    go_terms_json,
                    names_json,
                    accessions_json,
                    organisms_json,
                    length,
                    updated_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(uniref_id) DO UPDATE SET
                    go_terms_json = excluded.go_terms_json,
                    names_json = excluded.names_json,
                    accessions_json = excluded.accessions_json,
                    organisms_json = excluded.organisms_json,
                    length = excluded.length,
                    updated_at = excluded.updated_at
                ''',
                (
                    uniref_id,
                    json.dumps(go_terms),
                    json.dumps(names),
                    json.dumps(accessions),
                    json.dumps(organisms),
                    seq_len,
                    int(time.time()),
                ),
            )
            conn.commit()

    def delete_cache(self, uniref_ids: list[str] | None = None, delete_schema: bool = False) -> int:
        '''Delete cached UniRef entries.

        Parameters
        ----------
        uniref_ids : list of str or None, optional
            If None (default), delete all cached entries.
            If provided, delete only the matching UniRef IDs.
        delete_schema : bool, optional
            If True, also delete the cache database schema (i.e., drop the table). Default is False.

        Returns
        -------
        int
            Number of cache rows deleted.
        '''
        with sqlite3.connect(self.cache_db_path) as conn:
            if uniref_ids is None:
                cur = conn.execute('DELETE FROM uniref_cache')
                conn.execute('DELETE FROM go_term_cache')
            else:
                if not uniref_ids:
                    return 0
                placeholders = ','.join(['?'] * len(uniref_ids))
                cur = conn.execute(
                    f'DELETE FROM uniref_cache WHERE uniref_id IN ({placeholders})',
                    tuple(uniref_ids),
                )
            if delete_schema:
                if uniref_ids is not None:
                    raise ValueError('Cannot delete schema when specific uniref_ids are provided')
                conn.execute('DROP TABLE IF EXISTS uniref_cache')
                conn.execute('DROP TABLE IF EXISTS go_term_cache')
            conn.commit()
            return cur.rowcount

    def get_seq_annotation_strings(self, feature):
        '''Get nice string summaries of annotations for a given sequence

        Parameters
        ----------
        feature : str
            the uniref ID (i.e. 'UniRef90_A0A174LDE8') to query uniref about

        Returns
        -------
        shortdesc : list of (dict,str) (annotationdetails, annotationsummary)
            a list of:
                annotationdetails : dict containing at least the following:
                    'seqid' : str
                        the feature annotated (i.e. the sequence string)
                    'annotationtype : str
                        indicates the color for the list. we should improve this part in calour. currently set to "other"
                        and the list entry will be colored black
                    ADDITIONAL ITEMS:
                        for internal database interface use.
                        (for example details the database can use for the show_annotation_info() function called when an item in the list
                        is double clicked in the heatmap GUI)
                annotationsummary : str
                    a short summary of the annotation. Displayed in the heatmap gui in the listbox.
        '''
        logger.debug('getting uniref info for %s' % feature)
        shortdesc = []
        shortdesc.append( ({'unirefid':feature,'annotationtype':'other'},feature) )
        res = self._get_uniref_info(feature)
        if len(res['accessions']) == 0:
            logger.warning('uniref query failed:\n%s', res)
            shortdesc.append( ({'annotationtype':'other'},'uniref query failed. code %d' % res.code) )
            logger.debug(shortdesc)
            return shortdesc
        if 'length' in res and res['length'] is not None:
            shortdesc.append([{'annotationtype':'other'},'length: %d' % res['length']])
        if 'names' not in res:
            logger.debug('names empty')
            res['names'] = []
        for cres in res['names']:
            shortdesc.append([{'annotationtype':'other'},'name: %s' % cres])
        for cres in res['go_terms']:
            shortdesc.append([{'annotationtype':'other'},'go term: %s' % cres])
        for cres in res['organisms']:
            shortdesc.append([{'annotationtype':'other'},'organism: %s' % cres])
        logger.debug(shortdesc)
        return shortdesc

    def show_annotation_info(self, annotation):
        '''Show info about the annotation (can open an browser window or a python gui - whatever you prefer)
        Called when double clicking an annotation from the heatmap gui.
        annotation is the dict part returned from get_seq_annotation_strings().

        Parameters
        ----------
        annotation : dict
            should contain 'feature'
        '''
        # this is an example from the phenodb calour interface
        # should overwrite!
        uniref_base = 'https://www.uniprot.org/uniref/'
        webbrowser.open(uniref_base + annotation, new=2)


    def _get_quickgo_annotations(self, uniprot_accession: str, limit=200):
        """
        Fetch GO annotations from QuickGO for a UniProtKB accession.

        Parameters
        ----------
        uniprot_accession : str
            The UniProtKB accession (e.g., "P12345") to query for GO annotations.
        limit : int, optional
            The maximum number of GO annotations to retrieve (default is 200).

        Returns
        -------
        list of dict
            A list of GO annotations, where each annotation is represented as a dictionary containing:
                - "go_id": The GO term ID (e.g., "GO:0008150").
                - "aspect": The GO aspect (e.g., "biological_process", "molecular_function", "cellular_component").
                - "evidence_code": The evidence code supporting the annotation (e.g., "IDA", "IEA").
                - "go_name": The name of the GO term (e.g., "biological_process").
                - "assigned_by": The source that assigned the annotation (e.g., "UniProtKB").
        """
        if len(uniprot_accession)==0:
            return []
        url = f"{self.quickgo_url}/annotation/search"
        params = {
            "geneProductId": f"UniProtKB:{uniprot_accession}",
            "limit": limit,
            # 'includeFields': ["goId", "aspect", "evidenceCode", "goName", "assignedBy"]
        }
        logger.debug('querying quickgo for %s' % uniprot_accession)
        r = requests.get(url, params=params, headers={"Accept": "application/json"}, timeout=30)
        if r.status_code != 200:
            logger.warning('quickgo query (%s) failed with code %d. Error: %s' % (uniprot_accession, r.status_code, r.text))
            return []

        data = r.json()
        out = []
        for row in data.get("results", []):
            out.append({
                "go_id": row.get("goId"),
                "aspect": row.get("aspect"),
                "evidence_code": row.get("evidenceCode"),
                "go_name": row.get("goName"),
                "assigned_by": row.get("assignedBy"),
            })

        # deduplicate
        seen = set()
        dedup = []
        for x in out:
            key = (x["go_id"], x["evidence_code"], x["assigned_by"])
            if key not in seen:
                seen.add(key)
                dedup.append(x)
        return dedup

    def _get_go_term_names(self, go_ids: list[str]) -> dict[str, str | None]:
        '''Fetch GO term names from QuickGO for a list of GO IDs.

        Parameters
        ----------
        go_ids : list of str
            A list of GO term IDs (e.g., ["GO:0008150", "GO:0003674"]).
        
        Returns
        -------
        dict of str to str or None
            A mapping from GO term IDs to their names. If a GO term ID is not found, it will be mapped to None.
        '''
        if not go_ids:
            return {}

        unique_go_ids = sorted(set(go_ids))
        mapping = self._get_cached_go_term_names(unique_go_ids)

        missing_go_ids = [go_id for go_id in unique_go_ids if go_id not in mapping]
        if missing_go_ids:
            url = f"{self.quickgo_url}/ontology/go/terms/" + ",".join(missing_go_ids)
            r = requests.get(url, timeout=30)
            if r.status_code != 200:
                logger.warning('quickgo GO term query failed with code %d. Error: %s' % (r.status_code, r.text))
                fetched_mapping = {go_id: None for go_id in missing_go_ids}
            else:
                data = r.json()
                fetched_mapping = {}
                for term in data.get("results", []):
                    fetched_mapping[term["id"]] = term.get("name")

                # Keep misses as None so we do not re-query these IDs every time.
                for go_id in missing_go_ids:
                    fetched_mapping.setdefault(go_id, None)

            self._set_cached_go_term_names(fetched_mapping)
            mapping.update(fetched_mapping)

        return {go_id: mapping.get(go_id) for go_id in go_ids}
    
    def fill_cache(self, uniref_ids: list[str], force: bool = False) -> dict:
        '''Fill cache for a list of UniRef IDs.

        Parameters
        ----------
        uniref_ids : list of str
            UniRef IDs to cache.
        force : bool, optional
            If False (default), existing cached entries are reused.
            If True, entries are requeried and cache is updated.

        Returns
        -------
        dict
            Summary with number of processed IDs and failures.
        '''
        summary = {
            'requested': len(uniref_ids),
            'succeeded': 0,
            'failed': 0,
        }

        total = len(uniref_ids)
        logger.info('starting uniref cache fill: total=%d force=%s', total, force)
        log_every = 25

        for idx, uniref_id in enumerate(uniref_ids, start=1):
            try:
                res = self._get_uniref_info(uniref_id, use_cache=not force)
                if res.get('accessions'):
                    summary['succeeded'] += 1
                else:
                    summary['failed'] += 1
            except Exception:
                logger.exception('failed to fill cache for %s', uniref_id)
                summary['failed'] += 1

            if (idx % log_every == 0) or (idx == total):
                logger.info(
                    'uniref cache fill progress: %d/%d (succeeded=%d failed=%d)',
                    idx,
                    total,
                    summary['succeeded'],
                    summary['failed'],
                )

        logger.info(
            'finished uniref cache fill: requested=%d succeeded=%d failed=%d',
            summary['requested'],
            summary['succeeded'],
            summary['failed'],
        )

        return summary

    def _get_uniref_info(self, uniref_id: str, use_cache: bool = True) -> dict:
        '''Get go terms for a uniref id
        First queries uniref to get the list of member accessions, then queries each accession in quickgo to get go terms, then returns a dict of uniref_id -> list of go terms

        Parameters
        ----------
        uniref_id : str
            the uniref id to query (i.e. 'UniRef90_A0A174LDE8')

        Returns
        -------
        dict
            a dict with the following keys:
                'go_terms' : list of str
                    the list of go terms associated with this uniref id
                'accessions' : list of str
                    the list of accessions associated with this uniref id
                'names' : list of str
                    the list of protein names associated with this uniref id
                'length' : int
                    the average sequence length of the members of this uniref id
                'organisms' : list of str
                    the list of unique organism names associated with this uniref id
        '''
        uniref_id = uniref_id.split('.txt')[0] # in case it contains the .txt suffix from the feature name
        if use_cache:
            cached = self._get_cached_uniref_info(uniref_id)
            if cached is not None:
                logger.debug('cache hit for %s', uniref_id)
                return cached

        url = f"{self.uniref_url}/uniref/{uniref_id}/members/stream"
        r = requests.get(url, timeout=30)
        if r.status_code != 200:
            logger.error('uniref query (%s) failed with code %d', uniref_id, r.status_code)
            return {'go_terms': [], 'accessions': [], 'names': [], 'length': 0}
        r = r.json()['results']
        logger.debug('uniref query for %s returned %d members', uniref_id, len(r))
        seq_lengths = []
        seq_len = 0
        qnames = []
        protein_names = []
        organism_names = []
        qnames_for_go_terms = []
        for member in r:
            qname = None
            if 'proteinName' in member:
                protein_names.append(member['proteinName'])
            if 'accessions' in member and member['accessions']:
                qname = member['accessions'][0]
                qnames_for_go_terms.append(qname.split('_')[0])
            elif 'memberIdType' in member:
                if member['memberIdType'] == 'UniProtKB ID' or member['memberIdType'] == 'UniParc':
                    qname = member['memberId']
                if member['memberIdType'] == 'UniProtKB ID':
                    qnames_for_go_terms.append(member['memberId'].split('_')[0])
            if qname is None:
                    if 'uniparcid' in member:
                        qname = member['uniparcid']
            if qname:
                qnames.append(qname)
            if 'sequenceLength' in member:
                seq_lengths.append(member['sequenceLength'])
            if 'organismName' in member:
                organism_names.append(member['organismName'])
        if seq_lengths:
            seq_len = int(np.mean(seq_lengths))
        
        # get the go term ids for all qnames
        qnames_for_go_terms = list(set(qnames_for_go_terms))
        all_term_ids = self._get_quickgo_annotations(','.join(qnames_for_go_terms))

        # get all go term names using the cache to minimize queries to quickgo
        go_terms = []
        go_term_names = {}
        go_term_ids_to_search = []
        for go_annotation in all_term_ids:
            # get the go term id
            cterm = go_annotation.get('go_id', None)
            # if not cached - add to search list, otherwise add the name to the go_terms list
            # (we do this to minimize the number of queries to quickgo, which can be slow and rate limited. we query in batch for all missing terms at the end)
            if cterm not in go_term_names:
                go_term_ids_to_search.append(cterm)
            else:
                go_terms.append(go_term_names[cterm])
        # now search for all missing terms and update the cache and the go_terms list
        if len(go_term_ids_to_search)>0:
            go_term_ids_to_search = list(set(go_term_ids_to_search))
            new_terms = self._get_go_term_names(go_term_ids_to_search)
            go_term_names.update(new_terms)
            for cterm in go_term_ids_to_search:
                go_terms.append(go_term_names[cterm])
        # finally, remove duplicates from the go_terms list
        go_terms = list(set(go_terms))

        protein_names = list(set(protein_names))
        res = {}
        res['go_terms'] = go_terms
        res['accessions'] = qnames
        res['names'] = protein_names
        res['length'] = seq_len
        res['organisms'] = list(set(organism_names))
        self._set_cached_uniref_info(uniref_id, res)
        return res

    def get_organisms(self, exp: UniRefExperiment, inplace: bool = False) -> UniRefExperiment:
        '''Add "num_organisms" and "max_organism" columns to the feature metadata, where:
        - "num_organisms" is the number of unique organisms associated with each feature (uniref id)
        - "max_organism" is the organism that appears most frequently among the members of the uniref id
        
        Parameters
        ----------
        exp : calour.UniRefExperiment
            The UniRefExperiment containing the feature metadata to update.
        inplace : bool, optional
            If True, modify the feature metadata in place. If False, return a new calour Experiment with the updated metadata. Default is False.
        Returns
        -------
        calour.UniRefExperiment
            If inplace is False, returns a new UniRefExperiment with the updated feature metadata. If inplace is True, returns the same UniRefExperiment with modified feature metadata.
        '''
        if not inplace:
            exp = exp.copy()
        all_org_counts = []
        for crow in exp.feature_metadata.iterrows():
            org_counts = defaultdict(int)
            res = self._get_uniref_info(crow[0])
            corganisms = res.get('organisms', [])
            # go over corganisms and change some names:
            # - if starts with "Candidatus " - remove this prefix
            # if it is "human gut metagenome" - change to "metagenome"
            # if it starts with "uncultured " - remove this prefix
            # if it is "gut metagenome" - change to "metagenome"
            for i, org in enumerate(corganisms):
                if org.startswith('Candidatus '):
                    corganisms[i] = org[len('Candidatus '):]
                elif org == 'human gut metagenome':
                    corganisms[i] = 'metagenome'
                elif org.startswith('uncultured '):
                    corganisms[i] = org[len('uncultured '):]
                elif org == 'gut metagenome':
                    corganisms[i] = 'metagenome'

            row_org = [x.split(' ')[0] for x in corganisms]
            # count the number of times each organism appears in the row    org_counts = {}
            for org in row_org:
                org_counts[org] += 1
            all_org_counts.append(org_counts)
        org_num = [len(x) for x in all_org_counts]
        # get the organism with max count for each entry
        max_org = [max(x, key=x.get) if x else None for x in all_org_counts]
        exp.feature_metadata['num_organisms'] = org_num
        exp.feature_metadata['max_organism'] = max_org
        return exp
    
