# Appendix

## edge_path_to_paf_path
- [x] cut_loc 이름 바꾸기 -> edited_loc
- [x] Solve 함수의 경로 Return
    - [x] PafOutputData 설계
        - ctg_index
        - qry str, qry end
        - ref str, ref end
    - [x] Kpath 분석
        - 지금 원하는 경로에 대해서 PafOutputData들의 Sequence를 만들어주면 됨
        - K path를 호출하면 현재 반환하는건
        - std::vector<std::tuple<int, int, Distance>> path;
        - 각각 (a, b, w) 간선을 실제로 사용하고 있다는 것을 의미함
    - [x] `vector<PafOutputData> trace_path(const &path)` 이런식으로 각 path를 이제 recover해주는 함수를 설계할거임
        - path의 구성요소는 (a,b,w) 간선이므로
        - a,b의 형태에 따라서 Case work 를 해줘야함
        - 지금 총 3가지 형태가 있을 수 있음
        - 1. a가 start
            - 이 경우 항상 b는 (x, x) 형태이므로 그냥 path를 연결해주면 됨 ez
        - 2. a도 start가 아니고, b도 end가 아님
            - 2.1 a가 (x,x) 형태인 경우
                - b: (y,y)가 unoverlap이여서 (y,y) 형태면 (x,x) -> (y,y) 전이
                - b: (x,y)가 overlap이면 (x,x) -> (x,y) 전이
            - 2.2 a가 (x,y) 형태인 경우
                - b: (z,z) 또다른 overlap은 불가능하므로 (x,y) -> (z,z) 전이만이 가능함.
                - b: (y,z) 또다른 overlap을 만듦
        - 3. b가 end
            - 3.1 a: (x,x)
            - 3.1 a: (x,y)

- [x] K - Path Full Path에서 call_K_path를 true 상태 여러번 호출, 제대로 작동 X 문제
    - call_k_path일 경우 k+1개를 불러야 하는데 k를 불러서 생긴 문제였다.
    - 일단 default를 false로 바꾸고, k+1로 수정함.

## 구간 표현 방식
- 기본적으로 [st, en] 닫힌 구간 표현을 사용함

## Internal Vertex의 default vertex
- 기존 Internal Vertex는 처음에 Path를 찾을 때 쓰는 그래프 상의 노드.
- default vertex는 Internal Vertex를 default constructor로 만들 때만 true가 된다.
- default_vertex는 path를 다 찾고 나서, update_path를 할 때,
- default_vertex가 true면 ref 정보가 있을 때, false면 qry 정보만 있다.
