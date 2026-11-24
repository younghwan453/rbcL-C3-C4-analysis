# 연구 계획

## 질문과 사전 가설

| 연구 질문 | 사전 가설 또는 비교 기준 |
| --- | --- |
| C3와 C4에서 rbcL 전체 GC와 GC3 분포가 다른가? | 차이 방향은 사전에 지정하지 않습니다. GC3가 코돈 선택 차이에 더 민감할 수 있다는 가설을 검토합니다. |
| 코돈과 아미노산 조성에 그룹 차이가 있는가? | 같은 아미노산을 암호화하는 코돈은 여러 가지입니다. 코돈 조성 차이가 아미노산 변화 없이 관찰될 수 있습니다. |
| 문헌의 RbcL 149번과 309번 잔기가 표본에서도 광합성 유형과 연관되는가? | Flaveria에서 149 Ala와 309 Ile가 C4에 집중될 가능성을 검토합니다. 다른 계통에서는 패턴과 재현 여부를 별도로 확인합니다. |
| 관계가 외떡잎·진정쌍떡잎 및 가까운 분류군에서 어떻게 나타나는가? | 분류군별로 층화하고 속·과 내 사전 지정 비교를 제시합니다. |

## 입력 데이터와 수집 기준

1. 후보 CSV에 있는 정확한 학명으로 NCBI Taxonomy를 확인합니다. Taxonomy ID와 GenBank source의 학명 및 taxon 교차 참조를 대조합니다. 동의명은 원명과 대응 근거를 기록한 경우에만 허용합니다.
2. 각 그룹에서 한 종씩 시험합니다. C3 외떡잎은 Oryza sativa, C3 진정쌍떡잎은 Flaveria pringlei, C4 외떡잎은 Zea mays, C4 진정쌍떡잎은 Flaveria bidentis입니다. 네 종 통과 후 전체로 확대합니다.
3. 완전한 엽록체 RefSeq GenBank를 최우선으로 하고, 이어 완전 엽록체 GenBank와 주석이 충분한 full-length rbcL CDS를 검토합니다. 동일 우선순위의 복수 레코드는 종과 주석 일치 여부를 먼저 확인하고 선택 근거와 accession.version을 남깁니다. 종당 하나, accession은 중복 없이 유지합니다. 별개 종의 실제 동일 서열은 확인 후 각각 보존합니다.
4. CDS의 gene=rbcL, 위치, strand 및 완전성을 확인하고 feature 방향으로 추출합니다. 불확정 위치, partial, A/T/G/C 외 문자, 3의 배수가 아닌 길이, 유효하지 않은 시작·종결, 내부 stop 및 공개 번역과의 불일치를 검사합니다. CDS의 transl_table을 우선 적용하고 없으면 식물 색소체 table 11을 사용합니다. 주석 번역과 직접 번역 비교, protein accession, 검증 사유와 원본 GenBank를 남깁니다. [NCBI 유전 암호표](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11)
5. Biopython Entrez를 한 작업 흐름에서 이용하고 키 없이 초당 세 요청 이하로 유지합니다. 일시 오류는 제한적으로 재시도하고 실패·제외 사유를 기록합니다. NCBI_EMAIL은 환경변수에서 읽습니다. 자격 정보가 포함된 URL, 요청 내용은 출력하지 않습니다. [NCBI E-utilities 안내](https://www.ncbi.nlm.nih.gov/books/NBK25497/?report=printable)

## 주차별 처리와 출력

| 주차 | 처리 | 산출물 |
| --- | --- | --- |
| 12 | CDS·단백질 FASTA와 메타데이터를 동일 species_id로 연결. CDS 길이 및 염기 개수, AT·GC·GC1–3 계산. 광합성 유형과 분류군별 비교 | rbcL_nucleotide.fasta, rbcL_protein.fasta, species_metadata.csv, sequence_validation.csv, download_log.txt, 종별·그룹별 CSV와 GC 막대·상자·코돈 위치·산점도 PNG |
| 13 | 20개 표준 아미노산 개수와 백분율, ProteinAnalysis의 분자량·pI·GRAVY·aromaticity, 코돈 개수·빈도·RSCU | 조성·빈도·RSCU CSV, 아미노산·코돈·RSCU 히트맵, 그룹 막대·단백질 특성 상자 PNG |
| 14 | MUSCLE 단백질 정렬 후 코돈 역매핑. 보존·변이 위치, 쌍별 동일도, 잔기 빈도 및 후보 비교. CDS 정렬로 Neighbor-Joining | rbcL_aligned_protein.fasta, rbcL_aligned_cds.fasta, pairwise_identity.csv, amino_acid_substitutions.csv, candidate_residues.csv, phylogenetic_tree.nwk 및 C3/C4 색을 표시한 PNG |
| 15 | 종 단위 요약 및 탐색적 검정, 층화·가까운 분류군 비교, 다중검정 보정 | 표본 수·평균·표본 표준편차·효과크기·p-value 및 보정값을 포함한 표와 최종 보고서 |

전체 CDS의 길이와 A/T/G/C 원시 개수는 검증된 종결 코돈을 포함해 보고합니다. 비교용 AT·GC 백분율 및 GC1–3은 종결 코돈을 제거한 동일 CDS를 사용합니다. 전체 백분율 분모는 나머지 염기 수, 위치별 분모는 해당 위치의 염기 수로 각각 정의합니다. 코돈 빈도 분모는 유효 의미 코돈 수입니다. RSCU는 관측 코돈 수를 같은 아미노산의 동의 코돈에 균등 사용을 가정한 기대 개수로 나눕니다. 아미노산의 관측 수가 0이면 RSCU를 결측 처리합니다.

단백질 정렬에서 공백 열을 제외한 기준 서열의 잔기를 1부터 셉니다. Kapralov 등이 사용한 시금치(Spinacia oleracea) 구조 번호와 연결할 참조 RbcL을 14주차에 확인하고 출처, accession.version, 문헌 잔기 번호, 참조 단백질 잔기 번호, 정렬 열 및 종별 잔기를 기록합니다. 시금치는 번호 대응용이라 30종 집계 밖에 둡니다. 쌍별 동일도와 NJ의 비보정 p-distance는 둘 다 문자가 있는 위치를 분모로 하고 비교 가능 위치 수도 함께 기록합니다. gap이나 정보가 부족한 비교는 명시합니다. 잔기 간 차이와 분기 방향이 추정된 치환을 구분합니다.

Welch t-test는 그룹 평균 비교에 사용하고 효과크기는 평균차와 Hedges g를 제시합니다. 순위 비교가 적절한 경우 Mann–Whitney U와 순위 기반 효과크기를, 이분형 잔기는 Fisher 정확 검정과 오즈비를 제시합니다. 두 잔기 위치와 아미노산·코돈의 다중 비교는 각 사전 명시한 검정 묶음에 Benjamini–Hochberg FDR을 적용합니다. 네 그룹 할당은 설계로 정했으므로 이 표로 자연계 광합성 유형과 분류군의 연관을 검정하지 않습니다. 필요 시 조성 지표의 층화 차이와 기술적 유형×분류군 패턴을 기술합니다. 모든 종 단위 p-value는 계통학적 비독립성을 고려해 탐색적으로 해석하고, Flaveria·Atriplex·Cyperus와 Cleomaceae 내 비교를 함께 제시합니다.

## 오류 처리와 완료 기준

- 검증 실패는 레코드 단위로 남기고, 결측을 0으로 대체하거나 결과를 맞추기 위해 표본을 변경하지 않습니다. 교체가 필요하면 네 그룹 조건과 광합성 원문 근거를 다시 확인합니다.
- 정상·경계·오류·파일 연결을 포함한 테스트를 다섯 개 이상 실행하고, 손 계산용 예로 개수·GC·코돈 및 아미노산 결과를 대조합니다. 테스트 상태는 통과, 실패, 미실행으로 구분합니다.
- 각 주차에서 명령, 실제 실행 상태, 버전, 제외 이유, 출력과 한계를 기록합니다. 그래프는 CSV 입력과 수치가 일치하는지 확인하고 300 dpi PNG로 저장합니다.
- 최종 보고서는 Abstract, Introduction, Research Questions and Hypotheses, Materials and Methods, Results, Discussion, Limitations, Conclusion, References를 포함합니다. 깨끗한 작업 사본의 재현 실행에서 핵심 표와 경로를 확인합니다.
- 11주차 완료 조건은 문헌과 30종 후보 목록, 연구·학습 자료, 독립 환경과 버전 확인, 메타데이터 테스트 및 변경 검토입니다. 업로드는 별도 승인 후 진행합니다.
