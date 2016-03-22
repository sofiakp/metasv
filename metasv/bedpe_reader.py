import logging
import sys
import os

import vcf

from sv_interval import SVInterval

logger = logging.getLogger(__name__)

mydir = os.path.dirname(os.path.realpath(__file__))

bedpe_name = "BEDPE"
bedpe_source = set(["BEDPE"])

def extract_sv_info(filter_str, filter_names):
    filter_fields = filter_str.split(';')
    filter_vals = ['' for n in filter_names]
    for f in filter_fields:
        sub_filter_fields = f.split('=')
        if len(sub_filter_fields) < 2:
            continue
        for idx, name in enumerate(filter_names):
            if name == sub_filter_fields[0]:
                filter_vals[idx] = sub_filter_fields[1]
                if filter_vals[idx] == 'None':
                    filter_vals[idx] = None
    return filter_vals

class BedpeRecord:
    def __init__(self, record_string):
        self.name = bedpe
        fields = record_string.split()
        self.chr1 = fields[0]
        self.pos1 = int(fields[1])
        self.end1 = int(fields[2])
        self.chr2 = fields[3]
        self.pos2 = int(fields[4])
        self.end2 = int(fields[5])
        self.name = fields[6]
        sv_type = extract_sv_info(fields[11], ['TYPE'])[0]
        if sv_type == '':
            self.sv_type = 'UNK'
        else:
            self.sv_type = sv_type
        self.info = {
            "BP_CHR1": self.chr1,
            "BD_POS1": self.pos1,
            "BD_END1": self.end1,
            "BD_CHR2": self.chr2,
            "BD_POS2": self.pos2,
            "BD_END2": self.end2,
            "BD_NAME": self.name,
        }

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return "<" + self.__class__.__name__ + " " + str(self.__dict__) + ">"

    def to_sv_interval(self):
        if self.sv_type not in BedpeReader.svs_supported:
            return None

        if self.sv_type == "DEL" or self.sv_type == "INV":
            return SVInterval(self.chr1,
                              self.pos1 + 1,
                              self.end2
                              name=self.name,
                              sv_type=self.sv_type,
                              length=self.end2 - self.pos1,
                              sources=bedpe_pos,
                              cipos=[0, self.end1 - self.pos1],
                              info=self.info,
                              native_sv=self)
        else:
            logger.error("Bad SV type: " + repr(self))

        return None

    def to_vcf_record(self, sample):
        alt = [vcf.model._SV(self.sv_type)]
        sv_len = -self.sv_len if self.sv_type == "DEL" else self.sv_len
        info = {"SVLEN": sv_len,
                "SVTYPE": self.sv_type}
        if self.sv_type == "DEL" or self.sv_type == "INV":
            info["END"] = self.pos1 + self.sv_len
        else:
            return None

        info.update(self.info)

        vcf_record = vcf.model._Record(self.chr1,
                                       self.pos1,
                                       ".",
                                       "N",
                                       alt,
                                       ".",
                                       ".",
                                       info,
                                       "GT",
                                       [0],
                                       [vcf.model._Call(None, sample, vcf.model.make_calldata_tuple("GT")(GT="1/1"))])
        return vcf_record


class BedpeReader:
    svs_supported = set(["DEL"])

    def __init__(self, file_name, reference_handle=None, svs_to_report=None):
        logger.info("File is " + str(file_name))
        self.file_fd = open(file_name) if file_name is not None else sys.stdin
        self.reference_handle = reference_handle
        self.svs_supported = BedpeReader.svs_supported
        if svs_to_report is not None:
            self.svs_supported &= set(svs_to_report)

    def __iter__(self):
        return self

    def next(self):
        while True:
            line = self.file_fd.next().strip()
            if line:
                if line[0] != "#":
                    record = BedpeRecord(line)
                    if record.sv_type in self.svs_supported:
                        return record

    def get_header(self):
        return self.header
