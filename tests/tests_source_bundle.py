"""Offline scientific sources must survive a cache-free Git checkout."""
from pathlib import Path
import hashlib
import json
import tempfile
import unittest
from unittest.mock import patch
from atomic_database import SourceCache


class BundledSource(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name);self.bundle=self.root/'data';self.bundle.mkdir()
        self.url='https://example.org/pinned.dat';self.raw=b'1 2\n';self.name='pinned.dat'
        self.file=self.bundle/(hashlib.sha256(self.url.encode()).hexdigest()[:16]+'_'+self.name)
        self.file.write_bytes(self.raw);self.receipt=self.file.with_suffix('.dat.json')
        self.receipt.write_text(json.dumps(dict(url=self.url,sha256_bytes=hashlib.sha256(self.raw).hexdigest())),encoding='utf-8')

    def cache(self):return SourceCache(self.root/'runtime',offline=True,bundled_directory=self.bundle)

    def test_verified_bundle_is_read_without_network_or_copying_to_runtime(self):
        before={p.name:p.read_bytes() for p in self.bundle.iterdir()}
        with patch('atomic_database.urlopen',side_effect=AssertionError('network')):
            self.assertEqual(self.cache().get(self.name,dict(url=self.url,sha256=hashlib.sha256(self.raw).hexdigest())),self.raw)
        self.assertEqual(list((self.root/'runtime').iterdir()),[]);self.assertEqual(before,{p.name:p.read_bytes() for p in self.bundle.iterdir()})

    def test_corrupted_bundle_is_rejected(self):
        self.file.write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError,'Повреждён'):self.cache().get(self.name,dict(url=self.url))

    def test_bundle_requires_its_provenance_receipt(self):
        self.receipt.unlink()
        with self.assertRaisesRegex(ValueError,'квитанции'):self.cache().get(self.name,dict(url=self.url))

    def test_missing_sources_still_fail_offline(self):
        with self.assertRaisesRegex(OSError,'offline'):self.cache().get('missing',dict(url='https://example.org/missing'))
