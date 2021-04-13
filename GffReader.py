from os.path import abspath, expanduser
import gzip


class GffRecord(object):
    def __init__(self, seqid, source, type, start, end, score, strand, phase, attrs):
        try:
            assert "\n" not in seqid
            assert "\n" not in source
            assert "\n" not in type
            assert "\n" not in start
            assert "\n" not in end
            assert "\n" not in score
            assert "\n" not in strand
            assert "\n" not in phase
            assert "\n" not in attrs
        except AssertionError:
            raise ValueError("Invalid GFF record data")

        self.seqid = seqid
        self.source = source
        self._type = type
        self.start = int(start)
        self.end = int(end)
        self.length = self.end-self.start+1
        self.score = score
        self.strand = strand
        self.phase = phase
        self._attrs = attrs
        self.attributes = self._split_attr(attrs)

    def _split_attr(self, attributes):
        r = {}
        contents = attributes.split(";")
        for content in contents:
            if not content:
                continue
            if "=" not in content:
                print("%r is not a good formated attribute: no tag!") 
                continue
            tag, value = content.split("=", 1)
            r[tag] = value

        return r
    
    def to_string(self):
        attr = []
        for key, value in self.attributes.items():
            if key in "ID":
                attr.insert(0, "%s=%s" % (key, value))
            else:
                attr.append("%s=%s" % (key, value))
        r = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.seqid, self.source, self._type, self.start, self.end, self.score, self.strand, self.phase, ";".join(attr))
        return r

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, s):
        self._type = s    

    @classmethod
    def from_string(cls,s):
        try:
            assert "\n" not in s
            parts = s.split("\t")
            assert len(parts) == 9
            seqid, source, type, start, end, score, strand, phase, attributes = parts
           # assert strand in "+-"
        except AssertionError:
            raise ValueError("%r not recognized as a valid GFF record" % s)

        return GffRecord(seqid, source, type, start, end, score, strand, phase, attributes)


def open_gff(fn):
    filename = abspath(expanduser(fn))
    if filename.endswith(".gz"):
        ofs = gzip.open(filename, mode)
    elif filename.endswith(".dexta"):
        ofs = stream_stdout("undexta -vkU -w60 -i", filename)
    else:
        ofs = open(filename)

    for line in ofs.readlines():
        line = line.strip()
        if line.startswith('#'):
            continue
        if len(line) > 1:
            yield GffRecord.from_string(line)

    ofs.close()
    

