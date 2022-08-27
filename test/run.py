import os
import psycopg2
import sys
import traceback

_status = 0
_conn = psycopg2.connect(os.environ['DB_URI'])
_conn.cursor().execute("CREATE EXTENSION IF NOT EXISTS bioseqdb;")
_conn.commit()
def test(case):
    global _status
    try:
        case(_conn.cursor())
        print(f'{case.__name__} \x1B[1;32mok\x1B[0m', file=sys.stderr)
    except:
        traceback.print_exc()
        print(f'{case.__name__} \x1B[1;31merror\x1B[0m', file=sys.stderr)
        _status = 1
    _conn.rollback()

@test
def sequence_accept_basic(sql):
    sql.execute("SELECT 'ACGT'::NUCLSEQ;")
    assert sql.fetchone() == ('ACGT',)

@test
def sequence_accept_unknown(sql):
    sql.execute("SELECT 'N'::NUCLSEQ;")
    assert sql.fetchone() == ('N',)

@test
def sequence_accept_subsets(sql):
    sql.execute("SELECT 'RYKMSWBDHVN'::NUCLSEQ;")
    assert sql.fetchone() == ('RYKMSWBDHVN',)

@test
def sequence_reject_lowercase(sql):
    failed = False
    try:
        sql.execute("SELECT 'acgt'::NUCLSEQ;")
    except psycopg2.DataError as e:
        assert "invalid nucleotide in nuclseq_in: 'a'" in e.pgerror
        failed = True
    assert failed

_conn.close()
sys.exit(_status)
