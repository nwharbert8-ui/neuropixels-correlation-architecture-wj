"""
Download all 30 Neuropixels session NWB files from DANDI 000021.
Skips already-downloaded files. Downloads sequentially to avoid
overwhelming Drive sync.
"""
import os
import time
import urllib.request

DATA_DIR = r'G:\My Drive\inner_architecture_research\neuropixels_wj\data'
os.makedirs(DATA_DIR, exist_ok=True)

# All 30 sessions: (filename, asset_id, size_gb)
SESSIONS = [
    ('sub-707296975_ses-721123822.nwb', '224b57e5-c9a3-46ef-85db-966713f3ccbe', 1.6),
    ('sub-719828686_ses-754312389.nwb', '5a58bf3d-a1b9-444b-8ab0-ef5478aa42a6', 1.7),
    ('sub-740268983_ses-759883607.nwb', 'c9fc315f-5145-45a3-8183-3ad8fe499c4d', 1.7),
    ('sub-757329617_ses-773418906.nwb', '17a5432b-7193-42dd-a9ae-e9e8b3fef19c', 1.8),
    ('sub-739783158_ses-760345702.nwb', '47634abd-db85-48f5-9c33-01887a59d3bc', 1.8),
    ('sub-744915196_ses-762602078.nwb', 'cea32745-0d1a-4884-b06f-158d830bcbf3', 1.9),
    ('sub-722882751_ses-743475441.nwb', '522e8054-34ca-4579-ae80-350d0b24e0f4', 2.0),
    ('sub-716813540_ses-739448407.nwb', 'c7f32379-adce-4961-9211-d07790be9cab', 2.0),
    ('sub-718643564_ses-737581020.nwb', '106e84a7-1c41-43f1-a675-7df17e4aba69', 2.0),
    ('sub-735109599_ses-758798717.nwb', 'de872cf2-942d-4695-88b1-f7f89cb764fd', 2.0),
    ('sub-730760263_ses-755434585.nwb', 'edf10182-5a4c-454f-ad23-47987a5ca256', 2.1),
    ('sub-769360771_ses-791319847.nwb', '3adbba7c-3feb-468b-9829-33fcaa27aacd', 2.2),
    ('sub-730756767_ses-757970808.nwb', 'dfc3db15-066a-4a07-b615-a4d7e85c44e1', 2.2),
    ('sub-742602884_ses-763673393.nwb', '8a17b967-2aa9-4d6a-812c-92b62cf799d7', 2.4),
    ('sub-726170927_ses-746083955.nwb', '2e6df882-f31a-440b-a572-ba717a95bf80', 2.4),
    ('sub-734865729_ses-756029989.nwb', 'e0037b53-6d0e-4756-8ca4-f0dcb51ad339', 2.4),
    ('sub-719817799_ses-744228101.nwb', 'eb36f94f-d6e7-45c6-aa02-7d4ed23453d3', 2.5),
    ('sub-742714472_ses-761418226.nwb', 'cdeec90c-9ebe-431c-b3f7-b9f7d93ff37f', 2.5),
    ('sub-726298249_ses-754829445.nwb', '488a9bea-96fb-4028-a33a-08650f50be63', 2.5),
    ('sub-772616820_ses-799864342.nwb', '264f3b24-10d9-4100-98cf-d40b6c2a2068', 2.5),
    ('sub-726162193_ses-750749662.nwb', '3876c5f1-f38a-4c89-8f54-6128538f0066', 2.6),
    ('sub-699733573_ses-715093703.nwb', '58703c97-c0a9-4736-b684-73c85c1a444a', 2.7),
    ('sub-738651046_ses-760693773.nwb', '4513d0c9-1e2b-4c8a-aa22-ae81822537c9', 2.7),
    ('sub-745276220_ses-762120172.nwb', 'd6f3e82b-aca9-43fb-9810-048fc2124d50', 2.7),
    ('sub-717038285_ses-732592105.nwb', 'b4aeeb19-cdc6-4895-ab7b-bc8a688cf6f5', 2.7),
    ('sub-723627600_ses-742951821.nwb', '96c200cf-29c2-457a-b2f3-99f11de5b039', 2.7),
    ('sub-733457986_ses-757216464.nwb', '5f204c90-2005-4143-b880-07d80aa32b20', 2.8),
    ('sub-703279277_ses-719161530.nwb', '02291b99-e583-498b-9929-b68bba2c50e2', 2.9),
    ('sub-732548371_ses-751348571.nwb', 'be9f8fd8-8f16-4a66-acc6-9e04697650f3', 2.9),
    ('sub-726141242_ses-750332458.nwb', '286c7b06-3cde-4261-9090-e6fbe6c81945', 3.0),
]

print(f"Total sessions: {len(SESSIONS)}")
print(f"Total size: {sum(s[2] for s in SESSIONS):.1f} GB")

downloaded = 0
skipped = 0
failed = 0

for filename, asset_id, size_gb in SESSIONS:
    dest = os.path.join(DATA_DIR, filename)

    if os.path.exists(dest) and os.path.getsize(dest) > 1e9:
        actual_gb = os.path.getsize(dest) / (1024**3)
        print(f"SKIP {filename} (exists, {actual_gb:.1f} GB)")
        skipped += 1
        continue

    url = f'https://api.dandiarchive.org/api/assets/{asset_id}/download/'
    print(f"DOWNLOADING {filename} ({size_gb} GB)...", flush=True)
    t0 = time.time()

    try:
        urllib.request.urlretrieve(url, dest)
        elapsed = (time.time() - t0) / 60
        actual_gb = os.path.getsize(dest) / (1024**3)
        speed = actual_gb / (elapsed / 60) if elapsed > 0 else 0
        print(f"  Done in {elapsed:.1f} min ({actual_gb:.1f} GB, {speed:.1f} GB/hr)")
        downloaded += 1
    except Exception as e:
        print(f"  FAILED: {e}")
        failed += 1

print(f"\nSummary: {downloaded} downloaded, {skipped} skipped, {failed} failed")
print(f"Total on disk: {sum(os.path.getsize(os.path.join(DATA_DIR, f))/(1024**3) for f in os.listdir(DATA_DIR) if f.endswith('.nwb')):.1f} GB")
