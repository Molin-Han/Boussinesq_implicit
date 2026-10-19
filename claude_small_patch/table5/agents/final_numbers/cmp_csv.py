import json, csv
R = json.load(open('/home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/agents/final_numbers/myparse.json'))
rows = list(csv.DictReader(open('/home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/results.csv')))
seen = set()
bad = 0
for r in rows:
    t = r['tag']
    if t.startswith('old_'):
        print('old row', t, r['its'], r['pcapply_n'], r['pcapply_s'], r['np']); continue
    seen.add(t)
    if t not in R:
        print('csv tag not on disk', t); bad += 1; continue
    m = R[t]
    its_all = '/'.join(map(str, m['its']))
    pc_n = m['off_pc'][0] if m['off_pc'] else None
    pc_s = m['off_pc'][1] if m['off_pc'] else None
    ks = m['off_ksp'][1] if m['off_ksp'] else None
    patch = max(m['patch']) if m['patch'] else 0
    diffs = []
    if r['its_all'] != its_all: diffs.append(('its', r['its_all'], its_all))
    if str(pc_n if pc_n is not None else '') != r['pcapply_n']: diffs.append(('pcn', r['pcapply_n'], pc_n))
    if pc_s is not None and (r['pcapply_s'] == '' or abs(float(r['pcapply_s']) - pc_s) > 1e-9): diffs.append(('pcs', r['pcapply_s'], pc_s))
    if ks is not None and (r['ksp_s'] == '' or abs(float(r['ksp_s']) - ks) > 1e-9): diffs.append(('ksp', r['ksp_s'], ks))
    if str(patch) != r['patch'] and not (patch == 0 and r['patch'] in ('0', '')): diffs.append(('patch', r['patch'], patch))
    st = m['status'][:2] if m['status'] else []
    if st and (r['rc'] != st[0] or r['wall_s'] != st[1]): diffs.append(('status', r['rc'], r['wall_s'], st))
    if diffs:
        bad += 1; print('DIFF', t, diffs)
missing = sorted(set(R) - seen)
print('runs on disk not in csv:', missing)
print('n csv rows', len(rows), 'n disk', len(R), 'bad', bad)
