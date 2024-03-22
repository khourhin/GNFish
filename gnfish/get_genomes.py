#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import urllib
from Bio import Entrez
import re
import click
import pandas as pd
from loguru import logger
from pathlib import Path

import pdb


class GenomeDownloader(object):
    """Downloads genomes, RNA and proteins from NCBI."""

    # FIXME Does having the genomic option actually makes sense ?
    def __init__(self, working_dir, email, query, rna, protein, exclusive, retmax, refine):
        self.working_dir = Path(working_dir)
        params_to_log = {key: value for key, value in locals().items() if value is not None}
        logger.info(f"Running --{params_to_log} argument. Check help or README file for further information.")
        data_lst = []
        self.db = "assembly"
        self.data_dir = self.working_dir / "Data/"
        self.manage_create_directory()
        self.log_file = self.data_dir / "downloaded_genomes_log.tsv"
        self.log_file_data = self.check_file()
        self.query_lst = query.read().split("\n")
        self.email = email
        self.data_lst = ["genomic"]
        self.retmax = retmax
        self.refine = refine
        self.exclusive = exclusive
        if rna:
            data_lst.append("rna")
        if protein:
            data_lst.append("protein")

    # def create_directory(self, directory, path):
    #     directory = path + "/" + directory
    #     if not os.path.isdir(directory):
    #         os.mkdir(directory)
    #     else:
    #         logger.warning(f"{directory} directory already exits. New data will be added to previous one.")

    # FIXME path should not be returned anymore
    def manage_create_directory(self):
        self.data_dir.mkdir()
        (self.data_dir / "Genomic").mkdir()
        (self.data_dir / "Rna").mkdir()
        (self.data_dir / "Protein").mkdir()

    def check_file(self):
        if os.path.isfile(self.log_file):
            logger.warning(f"{self.log_file} file already exists.")
            logger.debug(
                f"New extraction results will be concatenated at the end of {self.log_file} and the already extant ones will be updated."
            )
            log_file_data = pd.read_csv(self.log_file, sep="\t")
            log_file_data["Assembly_ID"] = log_file_data["Assembly_ID"].astype(str)
        else:
            log_file_data = pd.DataFrame(columns=["Species", "Assembly_ID", "Genomic", "Rna", "Protein"])

        return log_file_data

    def get_record_entrez(self, final_query):
        Entrez.email = self.email
        handle = Entrez.esearch(db=self.db, term=final_query, retmax=self.retmax, sort="Significance")
        record = Entrez.read(handle)
        handle.close()
        return record

    def parse_ncbi_summary(self, id_num):
        esummary_handle = Entrez.esummary(db=self.db, id=id_num, report="full")
        summary = Entrez.read(esummary_handle, validate=False)

        species = summary["DocumentSummarySet"]["DocumentSummary"][0]["Organism"]
        species = re.sub(" ", "_", species)
        species = re.search("(.*?_.*?)_.*", species)
        species = species.group(1)

        url = summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]
        if not url:
            url = summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_GenBank"]
        if not url:
            logger.warning(f"No available {db} data for {species}. Check {species} in 'https://www.ncbi.nlm.nih.gov/'")
            return None, None
        else:
            return species, url

    def check_directory(self, directory):
        directory = directory
        directory_content = os.listdir(directory)
        if len(directory_content) == 0:
            os.rmdir(directory)

    def get_assembly_data(self, url, species, id_num, data_type, ext=".fna.gz"):
        label = os.path.basename(url)
        try:
            if data_type == "protein":
                ext = ".faa.gz"
                # FIXME
            link = os.path.join(url, label + "_" + data_type + ext)
            link = re.sub("\\\\", "/", link)
            directory = data_type.capitalize() + "/" + species
            dir_path = self.data_dir / directory
            dir_path.mkdir()
            urllib.request.urlretrieve(link, "%s/%s" % (dir_path, species + "_" + id_num + "_" + data_type + ext))
            return True
        except:
            dir_path = self.data_dir / directory
            self.check_directory(dir_path)
            logger.warning(
                f"No available {data_type} assembly data for {species} {id_num} or unable to download it. Check {species} {id_num} on 'https://www.ncbi.nlm.nih.gov/' assembly database."
            )
            return False

    def print_already_downloaded_info(self, data_typeo, species, id_num):
        logger.info(
            f"{data_type.capitalize()} data for {species} {id_num} is already downloaded. Check Data/{data_type.capitalize()} directory."
        )

    def print_successful_download(self, data_type, species, id_num):
        logger.info(
            f"{data_type.capitalize()} data for {species}_{id_num} was successfully downloaded. Check Data/{data_type.capitalize()} directory."
        )

    def launch_get_assembly_data(self, url, species, id_num):
        downloaded = False

        for data_type in self.data_lst:
            result = self.log_file_data.loc[self.log_file_data["Assembly_ID"] == id_num, data_type.capitalize()]
            if result.iloc[0] == 0:
                downloaded = self.get_assembly_data(url, species, id_num, data_type)
            elif result.iloc[0] > 0:
                self.print_already_downloaded_info(data_type, species, id_num)
            if (
                not self.exclusive
                and data_type != "genomic"
                and not downloaded
                and int(self.log_file_data.loc[self.log_file_data["Assembly_ID"] == id_num, "Genomic"].iloc[0]) == 0
            ):
                data_type = "genomic"
                downloaded = self.get_assembly_data(url, species, id_num, data_type)
            if downloaded:
                self.log_file_data.loc[self.log_file_data["Assembly_ID"] == id_num, data_type.capitalize()] = 1
                self.print_successful_download("genomic", species, id_num)

    def fetch_ncbi(self):
        """From parsed query file, fetch corresponding data from NCBI"""

        for query in self.query_lst:
            if query == "":
                continue
            final_query = f"{query} {self.refine}"
            logger.info(f"Searching for {final_query} at NCBI Assembly database.")
            record = self.get_record_entrez(final_query)
            # FIXME Change tyr/except to if statement
            try:
                error_sentence = record["ErrorList"]["PhraseNotFound"]
                logger.error(
                    f"{error_sentence} from {final_query} not found and can be unpredictable. Check out '{final_query}' on 'https://www.ncbi.nlm.nih.gov/' assembly database or correct {error_sentence} terms.\n"
                )
            except KeyError:
                ids = record["IdList"]
                for id_num in ids:
                    species, url = self.parse_ncbi_summary(id_num)
                    if species and url:
                        if id_num not in self.log_file_data["Assembly_ID"].values:
                            new_id_num = {
                                "Species": species,
                                "Assembly_ID": id_num,
                                "Genomic": 0,
                                "RNA": 0,
                                "Protein": 0,
                            }
                            self.log_file_data.loc[len(self.log_file_data)] = new_id_num
                            # FIXME Feasible with pathlib ?

                        self.launch_get_assembly_data(url, species, id_num)
                        self.log_file_data.to_csv(self.log_file, sep="\t", index=False)
                        logger.info(
                            "Genomes download ended. At Data/downloaded_genomes_log.tsv. You can find information about the downloaded genomes."
                        )


@click.command()
@click.argument("working_dir", type=click.Path(exists=True, dir_okay=True), required=True)
@click.argument(
    "email",
    type=click.STRING,
    required=True,
)
@click.argument(
    "query",
    type=click.File("r"),
    required=True,
)
@click.option("--rna", is_flag=True, help="Downloads rna annotation data.")
@click.option("--protein", is_flag=True, help="Downloads protein annotation data.")
@click.option(
    "--exclusive",
    is_flag=True,
    help="Download just protein or RNA annotation data if available.",
)
@click.option(
    "--retmax",
    type=click.INT,
    default=200,
    help="Number of NCBI records reported for every query.",
    show_default=True,
)
@click.option(
    "--refine",
    type=click.STRING,
    default="",
    show_default=True,
    is_flag=False,
    flag_value="AND (latest[filter] AND 'representative genome'[filter] AND all[filter] NOT anomalous[filter])",
    # FIXME Reformulate help
    help='Adds filter or field information to all queries. Flag value "AND (latest[filter] AND "representative genome"[filter] AND all[filter] NOT anomalous[filter])". Follow constant value structure for your custom refine.',
)
def main(working_dir, email, query, rna, protein, exclusive, retmax, refine):
    """Download genome, transcript, and protein data from NCBI database.

    WORKING_DIR set the path to you working directory where Data folder is going to be created.

            EMAIL your email for connecting to NCBI database.

        QUERY path to the file with your queries.
    """

    gd = GenomeDownloader(working_dir, email, query, rna, protein, exclusive, retmax, refine)
    gd.fetch_ncbi()
